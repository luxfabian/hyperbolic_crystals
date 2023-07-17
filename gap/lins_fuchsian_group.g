# -----------------------------------------------------------
#   ./gap/lins_fuchsian_group.g
#
#   Author: Fabian R. Lux
#   Date:   6/23/23
#
#   We perform a low-index normal subgroup search for the 
#   Fuchsian group of genus 2 as suggested in
#
#   Maciejko, Rajan, PNAS 119 (9) e2116869119 (2022)
#   https://doi.org/10.1073/pnas.2116869119
#
# -----------------------------------------------------------

Display("~~~~ GAP 4.12.1 ~~~~");

LoadPackage("LINS");

# -----------------------------------------------------------

Display("Constructing the {8,8} Fuchsian group:");

F := FreeGroup("g1","g2","g3","g4");

G := F / [ F.1 * F.2^(-1) * F.3 * F.4^(-1) * F.1^(-1) * F.2 * F.3^(-1) * F.4 ];

Display(G);

# -----------------------------------------------------------
n := 10;
Display("Generating normal subgroups up to index:");
Display(n);

# all normal subgroups up to n
lins_search := LowIndexNormalSubgroupsSearchForAll(G,n);
# lins_search := LowIndexNormalSubgroupsSearchForIndex(G,n,infinity);

# only normal subgroups of index equal to n
lins_nodes := ComputedNormalSubgroups(lins_search);

# converting LINS data structure to standard GAP data structure
subgroups := List( lins_nodes, node -> Grp(node)  );

# calculate the quotient groups
# quotients := List( List(subgroups), node -> StructureDescription(FactorGroup(G,node)));
# quotient_relations := List( subgroups, node -> FactorGroup(G,node));

# # calculate the group action on the quotient space
coset_tables := List( subgroups, node -> CosetTable(G,node));

PrintTo( "./gapout_fuchsian/coset_tables.csv" , coset_tables);
Display("Coset tables have been stored to ./gapout_fuchsian/coset_tables.csv");
# PrintTo( "./gapout/quotients.csv", quotients);
# Display("Classification of factor groups has been stored to ./gapout/quotients.csv");