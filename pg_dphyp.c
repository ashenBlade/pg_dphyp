#include "postgres.h"
#include "fmgr.h"
#include "nodes/bitmapset.h"
#include "optimizer/geqo.h"
#include "optimizer/paths.h"
#include "optimizer/pathnode.h"
#include "utils/builtins.h"
#include "utils/hsearch.h"
#include "utils/guc.h"
#include "miscadmin.h"

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

/* 
 * Represents a hypernode in query graph
 */
typedef struct HyperNode
{
	/* 
	 * Битовая карта отношений, которые этот узел представляет.
	 * Должен быть первым элементом в структуре.
	 */
	Bitmapset *nodes;

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

static join_search_hook_type prev_join_search_hook = NULL;

static bool enabled = true;

void _PG_init(void);
void _PG_fini(void);

static HyperNode *get_hypernode(DPHypContext *context, Bitmapset *nodes);

static inline int bms_first(Bitmapset *bms)
{
	int res = bms_next_member(bms, -1);
	Assert(res >= 0);
	return res;
}

static int us_leader(int *leaders, int a)
{
	if (leaders[a] == a)
		return a;
	else
		return leaders[a] = us_leader(leaders, leaders[a]);
}

static void us_union(int *leaders, int a, int b)
{
	int a_leader;
	int b_leader;

	a_leader = us_leader(leaders, a);
	b_leader = us_leader(leaders, b);

	leaders[a_leader] = b_leader;
}

static RelOptInfo *run_implicit_join_path(DPHypContext *context)
{
	/* 
	 * В запросе оказались неявные CROSS JOIN'ы (запятые без JOIN).
	 * Это может случиться и в случае подзапроса, где каждый JOIN предикат использует внешний параметр.
	 * В таких случаях, мы используем UNION-SET структуру для нахождения всех
	 * соединенных отношений (и получаем из RelOptInfo), а затем запускаем DPsize/GEQO для них.
	 */
	Bitmapset **disjoint_sets;
	List *disjoint_relations;
	ListCell *lc;
	int *leaders;
	int num_leaders;

	/* Инициализируем */
	num_leaders = list_length(context->initial_rels);
	leaders = palloc(sizeof(int) * num_leaders);
	for (size_t i = 0; i < num_leaders; i++)
	{
		leaders[i] = i;
	}
	
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
				Bitmapset *vertex = (Bitmapset *)lfirst(lc_vertex);
				int i = -1;
				while ((i = bms_next_member(vertex, i)) >= 0)
				{
					if (node->representative != i)
						us_union(leaders, node->representative, i);
				}
			}
		}
	}

	/* Для каждого лидера создаем его Bitmapset */
	disjoint_sets = palloc0(sizeof(Bitmapset *) * num_leaders);
	for (size_t i = 0; i < num_leaders; i++)
	{
		int leader = us_leader(leaders, i);
		if (disjoint_sets[leader] == NULL)
		{
			disjoint_sets[leader] = bms_make_singleton(i);
		}
		else
		{
			disjoint_sets[leader] = bms_add_member(disjoint_sets[leader], i);
		}
	}
	
	/* Заново проходимся и собираем все RelOptInfo */
	disjoint_relations = NIL;
	for (size_t i = 0; i < num_leaders; i++)
	{
		HyperNode *hypernode;

		if (disjoint_sets[i] == NULL)
			continue;
		
		hypernode = get_hypernode(context, disjoint_sets[i]);
		if (hypernode->rel == NULL)
			elog(ERROR, "failed to create RelOptInfo for disjoint set");
		
		disjoint_relations = lappend(disjoint_relations, hypernode->rel);
	}
	
	/* Запускаем обычный алгоритм JOIN для непересекающихся множеств */
	if (prev_join_search_hook)
		return prev_join_search_hook(context->root, list_length(disjoint_relations), disjoint_relations);
	else if (enable_geqo && list_length(disjoint_relations) >= geqo_threshold)
		return geqo(context->root, list_length(disjoint_relations), disjoint_relations);
	else
		return standard_join_search(context->root, list_length(disjoint_relations), disjoint_relations);
}

/* 
 * Создать Bitmapset со всеми битами от 0 до 'until' включительно
 */
static Bitmapset *create_all_bit_set(int until)
{
	Bitmapset *bv;

	bv = bms_make_singleton(until);
	for (int i = 0; i < until; ++i)
	{
		bv = bms_add_member(bv, i);
	}
	
	return bv;
}

static HTAB *create_dptable(List *base_hypernodes)
{
	ListCell *lc;
	HTAB *dptable;
	HASHCTL hctl;

	/*
	 * TODO: более умное определение размера таблицы - сейчас просто 1024 константа
	 */
	hctl.keysize = sizeof(Bitmapset *);
	hctl.entrysize = sizeof(HyperNode);
	hctl.hash = bitmap_hash;
	hctl.match = bitmap_match;
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
		entry->nodes = bms_copy(node->nodes);
	}

	return dptable;
}

static Bitmapset *get_neighbors(HyperNode *node, Bitmapset *excluded)
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
	Bitmapset *neighbors = NULL;

	/* Проходимся по всем гиперребрам */
	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;
		
		/* Проходимся по каждому узлу в гиперребре */
		foreach (lc2, edge)
		{
			Bitmapset *set = (Bitmapset *)lfirst(lc2);
			if (bms_overlap(set, node->nodes))
			{
				Bitmapset *me_excluded;
				Bitmapset *without_excluded;

				/* 
				 * Исключаю себя из Cross Join Set'а
				 */
				me_excluded = bms_difference(set, node->nodes);
				if (bms_is_empty(me_excluded))
					continue;

				/* 
				 * Если там кто-то остался, то нужно удалить исключенные
				 */
				without_excluded = bms_difference(me_excluded, excluded);
				if (bms_is_empty(without_excluded))
				{
					bms_free(me_excluded);
					continue;
				}

				/* Со всеми оставшимися узлами у меня есть JOIN ребро */
				neighbors = bms_union(neighbors, without_excluded);
				bms_free(me_excluded);
				bms_free(without_excluded);
			}
			else
			{
				int representative;

				if (bms_overlap(set, excluded))
					continue;
			
				Assert(!bms_is_empty(set));
				representative = bms_first(set);
				neighbors = bms_add_member(neighbors, representative);
			}
		}
	}
	
	return neighbors;
}

static bool hypernode_has_direct_edge_with(HyperNode *node, int id)
{
	ListCell *lc;

	Assert(!bms_is_member(id, node->nodes));

	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;
		foreach (lc2, edge)
		{
			Bitmapset *bms = (Bitmapset *)lfirst(lc2);
			
			/* Каждый узел - это cross join set, поэтому просто проверим, что этот узел есть в хоть одном узле ребра */
			if (bms_is_member(id, bms))
				return true;
		}
	}

	return false;
}

static bool hypernode_has_edge_with(HyperNode *node, Bitmapset *bms)
{
	ListCell *lc;

	Assert(!bms_overlap(node->nodes, bms));

	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;
		foreach (lc2, edge)
		{
			Bitmapset *cjs = (Bitmapset *)lfirst(lc2);
			if (bms_is_subset(cjs, bms))
				return true;
		}
	}

	return false;
}

static Bitmapset *create_bms(bitmapword *words, int nwords)
{
	Bitmapset *subset;

	subset = (Bitmapset *)palloc(offsetof(Bitmapset, words) + sizeof(bitmapword) * nwords);
	subset->type = T_Bitmapset;
	subset->nwords = nwords;
	
	for (size_t i = 0; i < nwords; i++)
	{
		subset->words[i] = words[i];
	}
	
	return subset;
}

static List *generate_all_subsets_long(Bitmapset *bms)
{
	/* 
	 * Эта функция используется, когда bms->words больше 1 слова, то есть алгоритм с простым вычитанием тут уже не сработает.
	 * В таком случае, мы разбиваем задачу по поиску комбинаций для 1 слова, на массив независимых слов.
	 * 
	 * То есть алгоритм следующий:
	 * 1. Храним массив bitmapset[]
	 * 2. Инициализируем его всеми 0
	 * 3. bitmapset[0] инициализируем начальным значением: (-init) & init
	 * 4. Начинаем алгоритм как для обычного перебора
	 * 5. Как только доходим до 0, то переходим к следующему слову
	 * 6. Для него запускаем 1 итерацию:
	 *      - Если 0, то (-init) & init (т.е. инициализируем)
	 *      - Иначе делаем 1 итерацию
	 * 7. Если после итерации получили 0, то переходим к следующему слову и делаем шаг 6 для него
	 * 8. Если на шаге 7 мы достигли последнего слова, то алгоритм заканчивается
	 * 
	 * Грубо говоря (очень), мы проводим выполняем сложение длинной арифметики:
	 * 1. Перебираем первое слово
	 * 2. Как только получили 0, то итерируем следующее слово (как бы перенос разряда после сложения)
	 * 3. Лавинно итерируем каждое очередное слово - переносим разряд, пока это возможно
	 * 4. Если добрались до последнего слова (для него надо сделать перенос) - конец
	 */
	bitmapword *words;
	List *subsets;
	int start_word;

	words = palloc0(sizeof(bitmapword) * bms->nwords);

	/* 
	 * Небольшая оптимизация - находим первое слово, которое гарантировано не 0
	 */
	start_word = 0;
	if (bms->words[0] == 0)
	{
		for (int i = 1; i < bms->nwords; i++)
		{
			if (bms->words[i] != 0)
			{
				start_word = i;
				break;
			}
		}
	}

	subsets = NIL;
	while (true)
	{
		int next_word;
	
		/* 
		 * Инициализируем первое слово: state = (-init) & init
		 */
		words[start_word] = (-bms->words[start_word]) & bms->words[start_word];

		/* 
		 * Итерируемся по нему до тех пор пока не получим 0 снова.
		 * Учитываем, что первое слово также может всегда быть 0.
		 */
		while (words[start_word] != 0)
		{
			subsets = lappend(subsets, create_bms(words, bms->nwords));

			words[start_word] = (words[start_word] - bms->words[start_word]) & bms->words[start_word];
		}

		/* 
		 * Теперь, words[start_word] == 0, поэтому переходим к следующему слову.
		 * При этом обрабатываем граничные случаи:
		 *  - начальное слово и есть последнее
		 *  - промежуточное слово может быть пустым
		 */

		if (start_word == (bms->nwords - 1))
			break;

		next_word = start_word + 1;
	
		while (next_word < bms->nwords)
		{
			while (next_word < bms->nwords && bms->words[next_word] == 0)
				++next_word;
			
			if (bms->nwords <= next_word)	
				break;

			if (words[next_word] == 0)
			{
				/* Если слово 0, то инициализируем его и начинаем новый цикл */
				words[next_word] = (-bms->words[next_word]) & bms->words[next_word];
				break;
			}

			/* В противном случае, делаем итерацию */
			words[next_word] = (words[next_word] - bms->words[next_word]) & bms->words[next_word];

			/* 
			 * В результате, у нас может получиться 0.
			 * Это значит, что текущее слово еще не до конца перебрали - продолжаем.
			 */
			if (words[next_word] != 0)
				break;

			/* Если текущее слово закончилось, то переходим к следующему */
			++next_word;
		}

		if (bms->nwords <= next_word)
			break;

		subsets = lappend(subsets, create_bms(words, bms->nwords));
	}

	pfree(words);
	return subsets;
}

static List *generate_all_subsets(Bitmapset *bms)
{
	bitmapword state;
	bitmapword init;
	List *result = NIL;

	/* 
	 * Используем тот алгоритм с генерацией битовой маски и постоянного вычитания: (prev - init) & init
	 */
	if (bms_is_empty(bms))
		return NIL;
	
	if (1 < bms->nwords)
		return generate_all_subsets_long(bms);

	Assert(bms->nwords == 1);

	init = bms->words[0];

	/* Получаем самый первый бит */
	state = (-init) & init;
	do
	{
		result = lappend(result, create_bms(&state, 1));
	} while ((state = (state - init) & init) != 0);

	return result;
}

static void emit_csg_cmp(DPHypContext *context, HyperNode *subgroup, HyperNode *complement)
{
	Bitmapset *nodes;
	RelOptInfo *joinrel;
	HyperNode *hypernode;

	/* 
	 * Возможно ли это? Может в таком случае просто возвращаться?
	 * Хотя, по построению алгоритма, у нас есть ребро, как минимум основание выполнить join
	 */
	Assert(subgroup->rel);
	Assert(complement->rel);

	nodes = bms_union(subgroup->nodes, complement->nodes);
	hypernode = get_hypernode(context, nodes);

	joinrel = make_join_rel(context->root, subgroup->rel, complement->rel);
	if (!joinrel)
		return;

	set_cheapest(joinrel);
	/* Лениво инициализируем */
	if (hypernode->rel == NULL)
        hypernode->rel = joinrel;
}

static void enumerate_cmp_recursive(DPHypContext *context, HyperNode *node, HyperNode *complement, Bitmapset *excluded)
{
	Bitmapset *complement_neighbors;

	List *neighbors_subsets;
	ListCell *lc;

	complement_neighbors = get_neighbors(complement, excluded);
	if (bms_is_empty(complement_neighbors))
		return;

	neighbors_subsets = generate_all_subsets(complement_neighbors);
	foreach(lc, neighbors_subsets)
	{
		Bitmapset *subset = (Bitmapset *)lfirst(lc);
		Bitmapset *neighbor_superset;
		HyperNode *superset_node;

		neighbor_superset = bms_union(complement->nodes, subset);
		superset_node = get_hypernode(context, neighbor_superset);

		if (superset_node->rel != NULL && 
			hypernode_has_edge_with(node, neighbor_superset))
		{
			emit_csg_cmp(context, node, superset_node);
		}

		bms_free(neighbor_superset);
	}
	/* 
	 * TODO: проследить где я создаю неправильный excluded
	 */

	/* 
	 * X = X u N(S_2, X)
	 */
	excluded = bms_union(excluded, complement_neighbors);
	list_free_deep(neighbors_subsets);

	bms_free(complement_neighbors);
	complement_neighbors = get_neighbors(complement, excluded);
	if (bms_is_empty(complement_neighbors))
	return;
	
	neighbors_subsets = generate_all_subsets(complement_neighbors);
	bms_free(complement_neighbors);

	foreach(lc, neighbors_subsets)
	{
		Bitmapset *subset = (Bitmapset *)lfirst(lc);
		Bitmapset *neighbor_superset;
		HyperNode *superset_node;

		neighbor_superset = bms_union(complement->nodes, subset);
		superset_node = get_hypernode(context, neighbor_superset);
		bms_free(neighbor_superset);

		enumerate_cmp_recursive(context, node, superset_node, excluded);
	}
}

static void emit_csg(DPHypContext *context, HyperNode *node)
{
	Bitmapset *excluded;
	Bitmapset *neighbors;
	int i;

	excluded = bms_union(node->nodes, create_all_bit_set(node->representative));
	neighbors = get_neighbors(node, excluded);
	if (bms_is_empty(neighbors))
	{
		bms_free(excluded);
		return;
	}

	i = -1;
	while ((i = bms_prev_member(neighbors, i)) >= 0)
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

static void enumerate_csg_recursive(DPHypContext *context, HyperNode *node, Bitmapset *excluded)
{
	Bitmapset *neighbors;
	List *subsets;
	Bitmapset *excluded_ext;
	ListCell *lc;

	neighbors = get_neighbors(node, excluded);
	if (bms_is_empty(neighbors))
		return;

	subsets = generate_all_subsets(neighbors);
	if (subsets == NIL)
		return;

	foreach(lc, subsets)
	{
		Bitmapset *neighbor_set = (Bitmapset *)lfirst(lc);
		Bitmapset *subset = bms_union(node->nodes, neighbor_set);
		HyperNode *subnode;
		
		subnode = get_hypernode(context, subset);
		if (subnode->rel != NULL)
		{
			Assert(subnode != NULL);
			emit_csg(context, subnode);
		}

		bms_free(subset);
	}

	/* X u N(S_1, X) */
	excluded_ext = bms_union(excluded, neighbors);

	foreach (lc, subsets)
	{
		Bitmapset *neighbor_set = (Bitmapset *)lfirst(lc);
		Bitmapset *subset = bms_difference(node->nodes, neighbor_set);
		HyperNode *subnode;

		subnode = get_hypernode(context, subset);
		enumerate_csg_recursive(context, subnode, excluded_ext);

		bms_free(subset);
	}

	bms_free(excluded_ext);
	list_free_deep(subsets);
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
		enumerate_csg_recursive(context, node, create_all_bit_set(i));
	}
}

static Bitmapset *map_to_internal_bms(List *initial_rels, Bitmapset *original)
{
	Bitmapset *target = NULL;
	ListCell *lc;
	int id = 0;

	foreach(lc, initial_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		if (bms_is_subset(rel->relids, original))
		{
			target = bms_add_member(target, id);
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
			Bitmapset *left_bms;
			Bitmapset *right_bms;

			left_bms = map_to_internal_bms(initial_rels, rinfo->left_relids);
			right_bms = map_to_internal_bms(initial_rels, rinfo->right_relids);
			if (!(bms_is_empty(left_bms) || bms_is_empty(right_bms)))
			{
				edge = list_make2(left_bms, right_bms);
			}
			else
			{
				bms_free(left_bms);
				bms_free(right_bms);
			}
		}
		else
		{
			/* Создаем CrossJoinSet */
			Bitmapset *bms;
			bms = map_to_internal_bms(initial_rels, rinfo->required_relids);
			if (!bms_is_empty(bms))
				edge = list_make1(bms);
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
				Bitmapset *bms;

				bms = map_to_internal_bms(initial_rels, em->em_relids);
				if (bms)
					edge = lappend(edge, bms);
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
	node->nodes = bms_make_singleton(id);

	return node;
}

static HyperNode *get_hypernode(DPHypContext *context, Bitmapset *nodes)
{
	HyperNode *node;
	Bitmapset *key = nodes;
	bool found;
	
	node = hash_search(context->dptable, &key, HASH_ENTER, &found);

	if (!found)
	{
		List *hyperedges;
		int id;

		node->nodes = bms_copy(nodes);
		node->representative = bms_first(nodes);
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
		while ((id = bms_next_member(nodes, id)) >= 0)
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
					Bitmapset *vertex = (Bitmapset *)lfirst(lc2);
					if (!bms_is_subset(vertex, nodes))
						new_edge = lappend(new_edge, vertex);
				}

				if (new_edge != NIL)
					hyperedges = lappend(hyperedges, new_edge);
			}
		}

		node->hyperedges = hyperedges;
	}

	return node;
}

static RelOptInfo *dphyp_join_search(PlannerInfo *root, int levels_needed, List *initial_rels)
{
	HTAB *dptable;
	HyperNode *result;
	Bitmapset *bms;
	DPHypContext context;
	List *base_hypernodes;
	ListCell *lc;
	int id;
	bool result_found;

	if (!enabled)
	{
		if (prev_join_search_hook)
			return prev_join_search_hook(root, levels_needed, initial_rels);
		if (enable_geqo && levels_needed >= geqo_threshold)
			return geqo(root, levels_needed, initial_rels);
		return standard_join_search(root, levels_needed, initial_rels);
	}

	/* 
	 * Шаги:
	 * 1. Каждому RelOptInfo назначаем ID (по факту индекс в массиве)
	 * 2. Создаем специальные структуры для представления гипер(узла)
	 * 3. Проходим joininfo и классы эквивалентности - создаем гиперребра (+ CrossJoinSet)
	 * 4. DPHyp
	 * 5. Пытаемся получить результат из всего BMS
	 * 6. Если есть - возвращаю
	 * 7. С помощью union-set нахожу все непересекающиеся множества
	 * 7. Для каждого получаю RelOptInfo из dptable
	 * 8. Запускаю DPsize старый (fallback на geqo)
	 */
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

	bms = create_all_bit_set(list_length(initial_rels) - 1);
	result = hash_search(dptable, &(bms), HASH_FIND, &result_found);
	if (result)
		return result->rel;
	
	/* 
	 * Union-set - находим все соединенные множества и получаем HyperNode для каждого из них
	 * затем запускаем старый алгоритм DPsize.
	 * 
	 * Это лучше делать с помощью DSU, но пока сделаю через Bitmapset - надо просто протестировать
	 */
	
	return run_implicit_join_path(&context);
}

void
_PG_init(void)
{

	DefineCustomBoolVariable("pg_dphyp.enabled",
							 "pg_dphyp join enumeration algorithm is enabled",
							 NULL,
							 &enabled,
							 enabled,
							 PGC_USERSET,
							 0, NULL, NULL, NULL);
	MarkGUCPrefixReserved("pg_dphyp");

	prev_join_search_hook = join_search_hook;
	join_search_hook = dphyp_join_search;
}

void
_PG_fini(void)
{
	join_search_hook = prev_join_search_hook;
}
