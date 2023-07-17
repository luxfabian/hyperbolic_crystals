# -----------------------------------------------------------
#   ./gap/lins_triangle_group.g
#
#   Author: Fabian R. Lux
#   Date:   6/23/23
#
#   We perform a low-index normal subgroup search for the 
#   proper triangle group similar to
#
#   Maciejko, Rajan, PNAS 119 (9) e2116869119 (2022)
#   https://doi.org/10.1073/pnas.2116869119
#
# -----------------------------------------------------------

Display("~~~~ GAP 4.12.1 ~~~~");

LoadPackage("LINS");

# -----------------------------------------------------------

Display("Constructing the proper {5,4} triangle group:");

F := FreeGroup("A","B");

p:=5;
q:=4;

G := F / [ F.1^p, F.2^q, (F.1*F.2)^2 ];

Display(G);

# -----------------------------------------------------------
n := 160;
Display("Generating normal subgroups up to index:");
Display(n);

# all normal subgroups up to n
lins_search := LowIndexNormalSubgroupsSearchForAll(G,n);
# lins_search := LowIndexNormalSubgroupsSearchForIndex(G,n,infinity);

# only normal subgroups of index equal to n
lins_nodes := ComputedNormalSubgroups(lins_search);

# converting LINS data structure to standard GAP data structure
subgroups := List( lins_nodes, node -> Grp(node)  );

generators := List( subgroups, node -> GeneratorsOfGroup(node) );

# group := subgroups[Length(subgroups)];

# newGenerators := GeneratorsOfGroup(group);

# newGenerators := RelatorsOfGroup(group);

Display(generators);

# calculate the quotient groups
quotients := List( List(subgroups), node -> StructureDescription(FactorGroup(G,node)));
quotient_relations := List( subgroups, node -> FactorGroup(G,node));

iso_quot := List( subgroups, node -> RelatorsOfFpGroup(SimplifiedFpGroup(Image( IsomorphismFpGroup(node) ))));

# iso_quot := List( quotient_relations, node -> SimplifiedFpGroup(node) );

# Display(iso_quot);

# # calculate the group action on the quotient space
coset_tables := List( subgroups, node -> CosetTable(G,node));

PrintTo( "./gapout_triangle/coset_tables.csv" , coset_tables);
Display("Coset tables have been stored to ./gapout_triangle/coset_tables.csv");
PrintTo( "./gapout_triangle/quotients.csv", quotients);
Display("Classification of factor groups has been stored to ./gapout_triangle/quotients.csv");