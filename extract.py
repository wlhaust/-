import random
from ete3 import Tree

t = Tree("RAxML_bipartitions_root.newick")

with open("result.csv", "a") as write_file:
    for each_node in t.traverse():
        if each_node.is_leaf():
            write_file.write(each_node.name.replace("'","") + "," + str(each_node.dist) + "\n")
