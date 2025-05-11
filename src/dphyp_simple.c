#include "postgres.h"
#include "utils/hsearch.h"
#include "optimizer/paths.h"
#include "optimizer/pathnode.h"
#include "miscadmin.h"

#include "dphyp_simple.h"
#include "unionset.h"
#include "simplebms.h"

/* 
 * Represents a hypernode in query graph
 */
typedef struct HyperNode
{
	/* 
	 * Битовая карта отношений, которые этот узел представляет.
	 * Должен быть первым элементом в структуре.
	 */
	bitmapword nodes;

	/* 
	 * Закешированный представитель 'nodes' (наименьший элемент)
	 */
	int representative;

	/* 
	 * Созданный RelOptInfo для него
	 */
	RelOptInfo *rel;

	/* 
	 * Список из гиперребер. Каждое гиперребро представляется как List из Bitmapset.
	 * Каждый Bitmapset - это гиперузел.
	 * При этом если длина внутреннего списка:
	 *
	 * 1 - CrossJoinSet
	 * 2 - обычное гиперребро
	 * >2 - мульти-гиперребро, по факту это multiway-join, но создан из классов эквивалентности
	 */
	List *hyperedges;
} HyperNode;

typedef struct DPHypContext
{
	/* 
	 * Информация о планировщике
	 */
	PlannerInfo *root;

	/* 
	 * Список из переданных RelOptInfo на вход
	 */
	List *initial_rels;

	/* 
	 * Таблица дин. прог.: Bitmapset -> HyperNode
	 */
	HTAB *dptable;

	/* 
	 * Список из изначальных гиперузлов.
	 */
	List *base_hypernodes;
} DPHypContext;

static HyperNode *get_hypernode(DPHypContext *context, bitmapword nodes);

static HTAB *create_dptable(List *base_hypernodes)
{
	ListCell *lc;
	HTAB *dptable;
	HASHCTL hctl;

	/*
	 * TODO: более умное определение размера таблицы - сейчас просто 1024 константа
	 */
	hctl.keysize = sizeof(bitmapword);
	hctl.entrysize = sizeof(HyperNode);
	hctl.hash = bmw_hash;
	hctl.match = bmw_match;
	hctl.hcxt = CurrentMemoryContext;
	dptable = (HTAB *)hash_create("DPHyp Hash Table", 1024L, &hctl,
								  HASH_ELEM | HASH_FUNCTION | HASH_COMPARE | HASH_CONTEXT);

	foreach (lc, base_hypernodes)
	{
		HyperNode *node = (HyperNode *)lfirst(lc);
		HyperNode *entry;
		bool found;

		entry = (HyperNode *) hash_search(dptable, &node->nodes, HASH_ENTER, &found);
		Assert(!found);

		entry->rel = node->rel;
		entry->hyperedges = node->hyperedges;
		entry->representative = node->representative;
		entry->nodes = node->nodes;
	}

	return dptable;
}

static bitmapword get_neighbors(HyperNode *node, bitmapword excluded)
{
	/* 
	 * Каждый узел в ребре обрабатываем независимо, в зависимости от того, есть ли там мой узел или нет.
	 *   Если меня в узле ребра нет, то:
	 *     1. Получаем представителя
	 *     2. Проверяем, что он не исключен
	 *   Если я есть, то:
	 *     1. Удаляю себя
	 *     2. Удаляю все исключенные
	 *     3. Добавляю все, что есть
	 */
	ListCell *lc;
	bitmapword neighbors;

	neighbors = 0;
	/* Проходимся по всем гиперребрам */
	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;
		
		/* Проходимся по каждому узлу в гиперребре */
		foreach (lc2, edge)
		{
			bitmapword set = (bitmapword )lfirst(lc2);
			if (bmw_overlap(set, node->nodes))
			{
				bitmapword me_excluded;
				bitmapword without_excluded;

				/* 
				 * Исключаю себя из Cross Join Set'а
				 */
				me_excluded = bmw_difference(set, node->nodes);
				if (bmw_is_empty(me_excluded))
					continue;

				/* 
				 * Если там кто-то остался, то нужно удалить исключенные
				 */
				without_excluded = bmw_difference(me_excluded, excluded);
				if (bmw_is_empty(without_excluded))
					continue;

				/* Со всеми оставшимися узлами у меня есть JOIN ребро */
				neighbors = bmw_union(neighbors, without_excluded);
			}
			else
			{
				int representative;

				if (bmw_overlap(set, excluded))
					continue;
			
				Assert(!bmw_is_empty(set));
				representative = bmw_first(set);
				neighbors = bmw_add_member(neighbors, representative);
			}
		}
	}
	
	return neighbors;
}

static bool hypernode_has_direct_edge_with(HyperNode *node, int id)
{
	ListCell *lc;

	Assert(!bmw_is_member(node->nodes, id));

	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;
		foreach (lc2, edge)
		{
			bitmapword bms = (bitmapword )lfirst(lc2);
			
			/* Каждый узел - это cross join set, поэтому просто проверим, что этот узел есть в хоть одном узле ребра */
			if (bmw_is_member(bms, id))
				return true;
		}
	}

	return false;
}

static bool hypernode_has_edge_with(HyperNode *node, bitmapword bms)
{
	ListCell *lc;

	Assert(!bmw_overlap(node->nodes, bms));

	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;
		foreach (lc2, edge)
		{
			bitmapword cjs = (bitmapword )lfirst(lc2);
			if (bmw_is_subset(cjs, bms))
				return true;
		}
	}

	return false;
}

static List *generate_all_subsets(bitmapword bms)
{
	bitmapword state;
	bitmapword init;
	List *result = NIL;

	/* 
	 * Используем тот алгоритм с генерацией битовой маски и постоянного вычитания: (prev - init) & init
	 */
	if (bmw_is_empty(bms))
		return NIL;


	init = bms;

	/* Получаем самый первый бит */
	state = (-init) & init;
	do
	{
		result = lappend(result, (void *) state);
	} while ((state = (state - init) & init) != 0);

	return result;
}

static void emit_csg_cmp(DPHypContext *context, HyperNode *subgroup, HyperNode *complement)
{
	bitmapword nodes;
	RelOptInfo *joinrel;
	HyperNode *hypernode;

	/* 
	 * Возможно ли это? Может в таком случае просто возвращаться?
	 * Хотя, по построению алгоритма, у нас есть ребро, как минимум основание выполнить join
	 */
	Assert(subgroup->rel);
	Assert(complement->rel);

	nodes = bmw_union(subgroup->nodes, complement->nodes);
	hypernode = get_hypernode(context, nodes);

	joinrel = make_join_rel(context->root, subgroup->rel, complement->rel);
	if (!joinrel)
		return;

	set_cheapest(joinrel);
	/* Лениво инициализируем */
	if (hypernode->rel == NULL)
        hypernode->rel = joinrel;
}

static void enumerate_cmp_recursive(DPHypContext *context, HyperNode *node, HyperNode *complement, bitmapword excluded)
{
	bitmapword complement_neighbors;

	List *neighbors_subsets;
	ListCell *lc;

	complement_neighbors = get_neighbors(complement, excluded);
	if (bmw_is_empty(complement_neighbors))
		return;

	neighbors_subsets = generate_all_subsets(complement_neighbors);
	foreach(lc, neighbors_subsets)
	{
		bitmapword subset = (bitmapword )lfirst(lc);
		bitmapword neighbor_superset;
		HyperNode *superset_node;

		neighbor_superset = bmw_union(complement->nodes, subset);
		superset_node = get_hypernode(context, neighbor_superset);

		if (superset_node->rel != NULL && 
			hypernode_has_edge_with(node, neighbor_superset))
		{
			emit_csg_cmp(context, node, superset_node);
		}
	}
	/* 
	 * TODO: проследить где я создаю неправильный excluded
	 */

	/* 
	 * X = X u N(S_2, X)
	 */
	excluded = bmw_union(excluded, complement_neighbors);
	list_free(neighbors_subsets);

	complement_neighbors = get_neighbors(complement, excluded);
	if (bmw_is_empty(complement_neighbors))
	return;
	
	neighbors_subsets = generate_all_subsets(complement_neighbors);

	foreach(lc, neighbors_subsets)
	{
		bitmapword subset = (bitmapword )lfirst(lc);
		bitmapword neighbor_superset;
		HyperNode *superset_node;

		neighbor_superset = bmw_union(complement->nodes, subset);
		superset_node = get_hypernode(context, neighbor_superset);

		enumerate_cmp_recursive(context, node, superset_node, excluded);
	}
}

static void emit_csg(DPHypContext *context, HyperNode *node)
{
	bitmapword excluded;
	bitmapword neighbors;
	int i;

	excluded = bmw_union(node->nodes, bmw_all_bit_set(node->representative));
	neighbors = get_neighbors(node, excluded);
	if (bmw_is_empty(neighbors))
		return;

	i = -1;
	while ((i = bmw_prev_member(neighbors, i)) >= 0)
	{
		HyperNode *complement;

		/* 
		 * В оригинальной статье здесь создается множество из единственного элемента.
		 * Затем в dptable ищется этот элемент и находится ребро, в которой правая часть - это подмножество множества выше.
		 * Но так как единственный возможный вариант - правая часть ребра - это и есть то самое множество: непустое подмножество множества из единственного элемента - это оно же.
		 */
		complement = (HyperNode *) list_nth(context->base_hypernodes, i);

		if (hypernode_has_direct_edge_with(node, i))
			emit_csg_cmp(context, node, complement);

		enumerate_cmp_recursive(context, node, complement, excluded);
	}
}

static void enumerate_csg_recursive(DPHypContext *context, HyperNode *node, bitmapword excluded)
{
	bitmapword neighbors;
	List *subsets;
	bitmapword excluded_ext;
	ListCell *lc;

	neighbors = get_neighbors(node, excluded);
	if (bmw_is_empty(neighbors))
		return;

	subsets = generate_all_subsets(neighbors);
	if (subsets == NIL)
		return;

	foreach(lc, subsets)
	{
		bitmapword neighbor_set = (bitmapword )lfirst(lc);
		bitmapword subset = bmw_union(node->nodes, neighbor_set);
		HyperNode *subnode;
		
		subnode = get_hypernode(context, subset);
		if (subnode->rel != NULL)
		{
			Assert(subnode != NULL);
			emit_csg(context, subnode);
		}
	}

	/* X u N(S_1, X) */
	excluded_ext = bmw_union(excluded, neighbors);

	foreach (lc, subsets)
	{
		bitmapword neighbor_set = (bitmapword )lfirst(lc);
		bitmapword subset = bmw_difference(node->nodes, neighbor_set);
		HyperNode *subnode;

		subnode = get_hypernode(context, subset);
		enumerate_csg_recursive(context, subnode, excluded_ext);
	}

	list_free(subsets);
}

static void solve(DPHypContext *context)
{
	int base_hypernodes_count = list_length(context->base_hypernodes);

	/* Оборачиваем список, чтобы итерироваться в обратном порядке */
	for (int i = base_hypernodes_count - 1; i >= 0; i--)
	{
		HyperNode *node = (HyperNode *) list_nth(context->base_hypernodes, i);
		CHECK_FOR_INTERRUPTS();
		emit_csg(context, node);
		enumerate_csg_recursive(context, node, bmw_all_bit_set(i));
	}
}

static bitmapword map_to_internal_bms(List *initial_rels, Bitmapset *original)
{
	bitmapword target;
	ListCell *lc;
	int id;

	target = 0;
	id = 0;
	foreach(lc, initial_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		if (bms_is_subset(rel->relids, original))
		{
			target = bmw_add_member(target, id);
		}

		++id;
	}
	
	return target;
}

static HyperNode *create_initial_hypernode(PlannerInfo *root, RelOptInfo *rel, int id, List *initial_rels)
{
	HyperNode *node;
	List *edges;
	ListCell *lc;

	node = (HyperNode *)palloc(sizeof(HyperNode));

	edges = NIL;

	foreach(lc, rel->joininfo)
	{
		RestrictInfo *rinfo = (RestrictInfo *)lfirst(lc);
		List *edge = NIL;

		/* 
		 * Такая проверка есть в make_restrictinfo_internal
		 */
		if (!bms_is_empty(rinfo->left_relids) &&
			!bms_is_empty(rinfo->right_relids) &&
			!bms_overlap(rinfo->left_relids, rinfo->right_relids))
		{
			/* Это бинарный оператор, поэтому создаем гиперребро */
			bitmapword left_bms;
			bitmapword right_bms;

			left_bms = map_to_internal_bms(initial_rels, rinfo->left_relids);
			if (bmw_is_empty(left_bms))
				continue;
			right_bms = map_to_internal_bms(initial_rels, rinfo->right_relids);
			if (bmw_is_empty(right_bms))
				continue;

			edge = list_make2((void *)left_bms, (void *)right_bms);
		}
		else
		{
			/* Создаем CrossJoinSet */
			bitmapword bms;
			bms = map_to_internal_bms(initial_rels, rinfo->required_relids);
			if (!bmw_is_empty(bms))
				edge = list_make1((void *)bms);
		}

		if (edge)
			edges = lappend(edges, edge);
	}

	/* 
	 * Из классов эквивалентности нельзя сделать CrossJoinSet,
	 * так как в каждом member может быть несколько соответствующих id
	 */
	if (rel->has_eclass_joins)
	{
		int i = -1;
		while ((i = bms_next_member(rel->eclass_indexes, i)) >= 0)
		{
			EquivalenceClass *ec = (EquivalenceClass *)list_nth(root->eq_classes, i);
			List *edge;

			/* 
			 * Ребра создаем только если в классе есть кто-то кроме самого себя
			 */
			if (bms_equal(ec->ec_relids, rel->relids))
				continue;

			edge = NIL;
			foreach(lc, ec->ec_members)
			{
				EquivalenceMember *em = (EquivalenceMember *)lfirst(lc);
				bitmapword bms;

				bms = map_to_internal_bms(initial_rels, em->em_relids);
				if (bms)
					edge = lappend(edge, (void *)bms);
			}

			if (edge != NIL)
			{
				edges = lappend(edges, edge);
			}
		}
	}

	node->hyperedges = edges;
	node->rel = rel;
	node->representative = id;
	node->nodes = bmw_make_singleton(id);

	return node;
}

static HyperNode *get_hypernode(DPHypContext *context, bitmapword nodes)
{
	HyperNode *node;
	bitmapword key = nodes;
	bool found;
	
	node = hash_search(context->dptable, &key, HASH_ENTER, &found);

	if (!found)
	{
		List *hyperedges;
		int id;

		node->nodes = nodes;
		node->representative = bmw_first(nodes);
		node->rel = NULL;

		hyperedges = NIL;

		/* 
		 * Создаем гиперребра для нового гиперузла.
		 * Проходим по всем базовым гиперузлам и каждому гиперребру.
		 * Если очередной узел этого ребра полностью находится в новом гиперузле,
		 * то просто не добавляем этот узел в ребро.
		 * 
		 * TODO: надо еще проверять на дубликаты
		 */
		id = -1;
		while ((id = bmw_next_member(nodes, id)) >= 0)
		{
			HyperNode *base_node = list_nth(context->base_hypernodes, id);
			ListCell *lc;

			foreach (lc, base_node->hyperedges)
			{
				List *edge = (List *)lfirst(lc);
				ListCell *lc2;
				List *new_edge = NIL;

				foreach(lc2, edge)
				{
					bitmapword vertex = (bitmapword )lfirst(lc2);
					if (!bmw_is_subset(vertex, nodes))
						new_edge = lappend(new_edge, (void *) vertex);
				}

				if (new_edge != NIL)
					hyperedges = lappend(hyperedges, new_edge);
			}
		}

		node->hyperedges = hyperedges;
	}

	return node;
}

static List *collect_disjoint_relations(DPHypContext *context)
{
	/* 
	 * В запросе оказались неявные CROSS JOIN'ы (запятые без JOIN).
	 * Это может случиться и в случае подзапроса, где каждый JOIN предикат использует внешний параметр.
	 * В таких случаях, мы используем UNION-SET структуру для нахождения всех
	 * соединенных отношений (и получаем из RelOptInfo), а затем запускаем DPsize/GEQO для них.
	 */
	List *disjoint_relations;
	List *disjoint_sets;
	ListCell *lc;
	unionset_state us_state;
	int num_leaders;

	/* Инициализируем */
	num_leaders = list_length(context->initial_rels);
	us_init(&us_state, num_leaders);
	
	/* Обхожу все гиперребра и составляю множества */
	foreach(lc, context->base_hypernodes)
	{
		HyperNode *node = (HyperNode *)lfirst(lc);
		ListCell *lc_edge;
		foreach(lc_edge, node->hyperedges)
		{
			ListCell *lc_vertex;
			List *edge = (List *)lfirst(lc_edge);

			foreach(lc_vertex, edge)
			{
				bitmapword vertex = (bitmapword )lfirst(lc_vertex);
				int i = -1;
				while ((i = bmw_next_member(vertex, i)) >= 0)
				{
					if (node->representative != i)
						us_union(&us_state, node->representative, i);
				}
			}
		}
	}

	/* Собираем все непересекающиеся множества */
	disjoint_sets = us_collect(&us_state);
	
	/* Заново проходимся и собираем все RelOptInfo */
	disjoint_relations = NIL;
	foreach(lc, disjoint_sets)
	{
		HyperNode *hypernode;
		Bitmapset *disjoint_set = (Bitmapset *)lfirst(lc);
		bitmapword simple_disjoint_set = 0;
		int i;

		i = -1;
		while ((i = bms_next_member(disjoint_set, i)) >= 0)
		{
			simple_disjoint_set = bmw_add_member(simple_disjoint_set, i);
		}
		
		hypernode = get_hypernode(context, simple_disjoint_set);
		if (hypernode->rel == NULL)
			elog(ERROR, "failed to create RelOptInfo for disjoint set");

		disjoint_relations = lappend(disjoint_relations, hypernode->rel);
	}

	list_free(disjoint_sets);
	us_free(&us_state);

	return disjoint_relations;
}


List *dphyp_simple(PlannerInfo *root, int levels_needed, List *initial_rels)
{
	HTAB *dptable;
	HyperNode *result;
	bitmapword all_rels_bit_set;
	DPHypContext context;
	List *base_hypernodes;
	ListCell *lc;
	int id;
	bool result_found;

	base_hypernodes = NIL;
	id = 0;
	foreach(lc, initial_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		HyperNode *node;

		node = create_initial_hypernode(root, rel, id, initial_rels);
		base_hypernodes = lappend(base_hypernodes, node);
		++id;
	}

	dptable = create_dptable(base_hypernodes);
	
	/* 
	 * Access paths for all rels in 'initial_rels' already found.
	 * Ready to run DPHyp.
	 */
	context.dptable = dptable;
	context.initial_rels = initial_rels;
	context.root = root;
	context.base_hypernodes = base_hypernodes;
	solve(&context);

	all_rels_bit_set = bmw_all_bit_set(list_length(initial_rels) - 1);
	result = hash_search(dptable, &(all_rels_bit_set), HASH_FIND, &result_found);
	if (result && result->rel)
		return list_make1(result->rel);

	/* 
	 * Union-set - находим все соединенные множества и получаем HyperNode для каждого из них
	 * затем запускаем старый алгоритм DPsize.
	 * 
	 * Это лучше делать с помощью DSU, но пока сделаю через Bitmapset - надо просто протестировать
	 */
	
	return collect_disjoint_relations(&context);
}
