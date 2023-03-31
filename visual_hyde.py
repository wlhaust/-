import os 
import sys
import argparse
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, faces, random_color, TextFace
from PIL import Image
matplotlib.use('Agg')
os.environ ['QT_QPA_PLATFORM'] ='offscreen'

Description = (
  '''
   ----------------------------------------------------------------------------- 
  |                              visual_hyde.py                                |
   ----------------------------------------------------------------------------- 

  Created by Jian He (j.he930724@gamail.com)                      
  This is version 1.1 created on 2022.4.21                                                                                                                            
                                                                                        
  Dependancies: 
  python3
  ete3 (conda install -c etetoolkit ete3 ete_toolchain)
  matplotlib (conda install matplotlib)
  numpy (conda install numpy)
  pandas (conda install pandas)

  Description:
  HyDe (Blischak et al., 2018) is a python package to detect hybridization using 
  phylogenetic invariants based on the coalescent model. Similar to Patterson's 
  D-statistic (Patterson et al. 2012), HyDe considers a rooted, 4-taxon tree 
  “((P1, H, P2), O)” that consisting of an outgroup “O” and a triplet of 
  ingroup. The hybrid taxa “H” is modeled as a hybrid species between two 
  parental taxa, “P1” and “P2”. According to Meng and Davey's hybrid model, the 
  hybrid taxa is either sister to “P1” with probability (γ) or sister to “P2” 
  with probability (1–γ), The null hypothesis was that when hybrid was absent, 
  the expected value of γ should be 0. HyDe perform a formal statistical test 
  of γ = 0 versus γ > 0 using Z-test. The higher the Z-score, the more reliable 
  of a hybrid event.
  
  To visualize the HyDe results, we draw the heatmap for each potential hybrid 
  sample “H”. For each figure, the parents (“P1” and “P2”) were placed in the 
  heatmap on “x” or “y” axes respectively. all the taxa on the axes were 
  arranged phylogenetically. Each square in the heatmap was represents a HyDe 
  test which preform in a 4-taxon tree “((P1, H, P2), O)”. The different color 
  of the squares corresponds different "γ" values that represent genetic 
  contribution of taxon on the left axes. 

  However, when a signal shows that the "H" was detected as a hybrid of "P1" 
  and "P2", it does not only represent a recent hybridization event that related 
  to the three taxa, but also possibly related with their respective ancestors. 
  Specifically, if the ancestor of "P1" or "P2" was involved in a hybridization 
  event, then all the descendants of the ancestor should be able to detect the 
  hybridization event also. Therefore, in a heatmap, if some squares with 
  similar colors formed a larger one in the heatmap, and their corresponding 
  taxa on the axis could form a monophyletic clade, then the most recent common 
  ancestor (MRCA) of the monophyletic clade may be indicated as the parent 
  of "H". 

  Further, if the MRCA of a clade was origin from a hybridization event, then 
  all the taxa on this clade should also detect the same hybridization event 
  signal. Therefore, for each clade, we stacked all the heatmaps which draw for 
  the leaves on this clade. For the heatmap after stacking, only if all the 
  stacked squares were light, the corresponding square will be light. In other 
  words, the light squares means that the hybridization event occurred on the 
  MRCA of these samples which involved in stacking.
  ''')


#如果用户没有给预先定义的大分枝，则简单定义一下大分枝，每个大分枝内的物种数量不超过5个
def make_predefined_clade_file(tree_file):
  def return_species_name_clade(t):
    child_node = t.children
    for each_node in child_node:
      if len(each_node.get_leaf_names()) <= maxnum_in_predefined_clade:
        
        for each_species_name in each_node.get_leaf_names():
          write_file.write(each_species_name + ",")
        write_file.write("\n")
      else:
        return_species_name_clade(each_node)
  t = Tree(tree_file)
  maxnum_in_predefined_clade = 5
  with open("Predefined_clade.txt", "w") as write_file:
    return_species_name_clade(t)

#处理树，得到一些信息
def parse_tree(tree_file, Predefined_clade_file): 
  #读取树。
  t = Tree(tree_file)
  Label = t.get_leaf_names()
  #获取排列好的分枝
  subtrees = []
  with open(Predefined_clade_file, "r") as read_file:
    for each_line in read_file:
      line_list = []
      for each_name in each_line.split(","):
        if len(each_name) > 2:
          line_list.append(each_name.replace("\n","").replace(" ",""))
      if line_list:
        line_list.append(Label.index(line_list[0]))
        subtrees.append(line_list)
    #按照Label中出现的先后顺序排序    
    def takeLast(elem):
      return elem[-1]
    subtrees.sort(key = takeLast)
    Highlight_subtrees = []
    for each_list in subtrees:
      Highlight_subtrees.append(each_list[:-1])
  #输出名称长度的列表
  name_len_list = []
  for each_name in Label:
    name_len_list.append(len(each_name))
  name_len_list.sort()
  name_len = (name_len_list[0], name_len_list[-1])
  return t, Label, Highlight_subtrees, name_len

#建立某个sample的hyde表格，表格的横纵坐标轴为sample name，表格的值为gamma
def make_hotmap_table_gamma(Label, hypothesis_hybrid_species, csv_file_name, 
                            zscore):
  #该函数通过Hyde结果获取三个种杂交验证可信度的γ值
  def get_gamma(hyde_out_table, each_index, each_column):
    sub_table = hyde_out_table[(hyde_out_table["P1"] == each_index) & 
                               (hyde_out_table["P2"] == each_column)]
    if list(sub_table["Gamma"]) == []:
      sub_table2 = hyde_out_table[(hyde_out_table["P1"] == each_column) & 
                                  (hyde_out_table["P2"] == each_index)]
      if list(sub_table2["Gamma"]) == []:
        pass
      else:
        return 1 - list(sub_table2["Gamma"])[0]
    else:
      return list(sub_table["Gamma"])[0]


  #先生成行和列数目是物种数目但内容随机的pandas表格，之后再将这些表格获取为相对应的γ值
  hyde_table = pd.read_csv(csv_file_name, sep="\t")
  sub_table = hyde_table[(hyde_table["Hybrid"] == hypothesis_hybrid_species) & 
                         (hyde_table["Zscore"] > zscore)] #筛选z值
  df = pd.DataFrame(np.random.randn(len(Label), len(Label)), index=Label, 
                    columns=Label)
  
  for each_index in Label:
    for each_column in Label:
      gamma = get_gamma(sub_table, each_index, each_column)
      if gamma != None:
        '''
        if gamma < 0.5:
          gamma = gamma - 0
        else:
          gamma = 1 - gamma
        '''
        df.at[each_index, each_column] = gamma
      else:
        df.at[each_index, each_column] = 0
  #下面的代码使热图只显示半个三角
  len_Label = len(Label)
  for each_num_index in range(len_Label):
    for each_num_index2 in range(len_Label):
      if each_num_index2 > each_num_index:
        df.iloc[each_num_index:each_num_index + 1, 
                each_num_index2:each_num_index2 + 1] = None  
  
  return df

#将所有待检测物种的热图表格都加入一个字典中
def make_hyde_output_array_dict(Label, csv_file_name, zscore):
  hyde_output_array_dict = {}
  n = 0
  for each_hypothesis_hybrid_species in Label:
    n = n + 1
    print("start parsing " + each_hypothesis_hybrid_species + 
          "'s hyde output " + str(n) + "/" + str(len(Label)))
    hyde_output_array = make_hotmap_table_gamma(Label, 
                          each_hypothesis_hybrid_species, csv_file_name, zscore)
    hyde_output_array_dict[each_hypothesis_hybrid_species] = hyde_output_array
  return hyde_output_array_dict

#找到某特定节点之后所有物种所共享的杂交事件，并输出一个表格
def find_common_hybrid_in_nodes(hyde_output_array_dict, each_node, Label):
  def get_num(each_index, each_column):
    gamma_list = []
    base_species_have_zscore = True

    gamma_left_child_list = []
    for each in each_node.children[0].get_leaf_names():
      gamma_left_child = hyde_output_array_dict[each].at[each_index, 
                                                         each_column]
      gamma_left_child_list.append(gamma_left_child)
    if sum(gamma_left_child_list) == 0:
      base_species_have_zscore = False
    gamma_right_child_list = []
    for each in each_node.children[1].get_leaf_names():
      gamma_right_child = hyde_output_array_dict[each].at[each_index, 
                                                          each_column]
      gamma_right_child_list.append(gamma_right_child)
    if sum(gamma_right_child_list) == 0:   
      base_species_have_zscore = False

    number_of_non_gamma_species = 0
    for each_leaves in each_node.get_leaf_names():
      gamma = hyde_output_array_dict[each_leaves].at[each_index, each_column]
      if gamma != 0:
        gamma_list.append(gamma)
      else:
        number_of_non_gamma_species = number_of_non_gamma_species + 1
    if base_species_have_zscore:
      if number_of_non_gamma_species/len(each_node.get_leaf_names()) <= 0.5:
        return sum(gamma_list)/(len(gamma_list) + 0.000001)
      else:
        return 0
    else:
      return 0
  df = pd.DataFrame(np.random.randn(len(Label), len(Label)), index=Label, 
                    columns=Label)
  for each_index in Label:
    if each_index not in each_node.get_leaf_names():
      for each_column in Label:
        if each_column not in each_node.get_leaf_names():
          gamma = get_num(each_index, each_column)
          df.at[each_index, each_column] = gamma
        else:
          df.at[each_index, each_column] = 0
    else:
      for each_column in Label:
        df.at[each_index, each_column] = 0
  #下面的代码使热图只显示半个三角
  len_Label = len(Label)
  for each_num_index in range(len_Label):
    for each_num_index2 in range(len_Label):
      if each_num_index2 > each_num_index:
        df.iloc[each_num_index:each_num_index + 1, 
                each_num_index2:each_num_index2 + 1] = None  
  return df  

#绘制热图
def draw_hotmap(pic_name, hyde_output_array):
  hyde_output_array.to_csv(pic_name + ".csv",index =True ,sep = ',')
  #创建空背景
  fig = plt.figure(figsize=(30,30))

  #设定axes大小，并将其填入fig中
  border_width = 0.00001
  ax_size = [0+border_width, 0+border_width, 
              1-2*border_width, 1-2*border_width]  
  ax = fig.add_axes(ax_size)

  #画出热图
  cmap = []
  color = 1
  lucency = 0
  for each in range(5000):
    cmap.insert(each, np.array([0,0,color,lucency]))
    color = color - (1/5000)
    lucency = lucency + (1/5000)
  color = 0
  lucency = 1
  for each in range(5001, 9999):
    cmap.insert(each, np.array([color,0,0,lucency]))
    color = color + (1/5000)
    lucency = lucency - (1/5000)
  newcmp = ListedColormap(cmap)
  im = ax.imshow(hyde_output_array, norm = 
                 matplotlib.colors.Normalize(vmin=0, vmax=1), cmap = newcmp)
  position=fig.add_axes([0.9, 0.2, 0.05, 0.7])
  plt.colorbar(im, cax=position)
  cbar = plt.colorbar(im, cax=position)
  #labels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
  cbar.ax.tick_params(labelsize=50)
  #norm = matplotlib.colors.Normalize(vmin=0, vmax=0.5) 
  #热图的小方格间加入小空隙
  ax.set_xticks(np.arange(hyde_output_array.shape[1]+1)-.5, minor=True)
  ax.set_yticks(np.arange(hyde_output_array.shape[0]+1)-.5, minor=True)
  ax.grid(which="minor", color="black", linestyle='-', linewidth=2)

  #保存热图至"hotmap.png"中
  plt.savefig("hotmap.png", dpi = 200)

  #清理掉内存中的图
  plt.cla()
  plt.close("all")

#绘制树
def draw_tree(tree_file, node_num, highlight_clade, Clade_file, name_len, 
              picture_size):

  #重新获取物种树
  t = Tree(tree_file)

  #该函数可绘制不需要高亮显示的内部分枝
  def node_layout_background(node):
    ns = NodeStyle()
    ns["hz_line_width"] = 4.5
    ns["vt_line_width"] = 4.5
    ns["size"] = 0
    node.set_style(ns)  

  #该函数可绘制需要高亮显示的大分枝
  def node_layout(node, color):
    if node.is_leaf():
      descFace = faces.TextFace(node.name + 
                                "·"*((name_len[1]-len(node.name)) * 2 + 10), 
                                fsize=30, fgcolor = color)
      descFace.margin_top = 3
      descFace.margin_bottom = 3
      descFace.border.margin = 1
      node.add_face(descFace, column=1, position='aligned') 
    ns = NodeStyle()
    ns["hz_line_width"] = 4.5
    ns["vt_line_width"] = 4.5
    ns["vt_line_color"] = color
    ns["hz_line_color"] = color  
    ns["size"] = 0
    node.set_style(ns)

  #首先绘制内部不许要高亮的节点以及不在用户所定义的分枝列表中的节点。
  Clade_file_temp_list = [] #先获取用户所定义的分枝列表中的所有叶子
  for each1 in Clade_file:
    for each2 in each1:
      Clade_file_temp_list.append(each2)
  for each_node in t.traverse():
    #这里的逻辑如下：遍历物种树内部的所有节点，如果是内部节点，则设定树分支的格式。
    #如果一个节点是叶子节点，如果这个节点是正在检测的物种，则高亮该分支，
    #如果不在用户所定义的大分枝中，则把他画成黑色
    if each_node.is_leaf():
      if each_node.get_leaf_names()[0] == highlight_clade:
        node_face = TextFace("o",fsize = 40)
        node_face.background.color = "red"
        each_node.add_face(node_face, column=0, position = "branch-right")
      if each_node.get_leaf_names()[0] not in Clade_file_temp_list:
        color = None
        node_layout(each_node, color)
    else:
      node_layout_background(each_node)

  #之后绘制用户所定义的大分枝。
  h = 0
  for each_subtree in Clade_file:
    color = random_color(h, s = 0.9, l = 0.4)
    h = h + 0.58
    if len(each_subtree) != 1:
      Highlight_node = t.get_common_ancestor(each_subtree)
      for each_node in Highlight_node.traverse():
        node_layout(each_node, color)
    else:
      Highlight_node = t.get_leaves_by_name(each_subtree[0])[0]
      node_layout(Highlight_node, color)

  #定义一些其他的树形
  node_face = TextFace(str(node_num),fsize = 40)
  node_face.background.color = "LightGreen"
  t.get_common_ancestor(highlight_clade).add_face(node_face, column=0, 
                        position = "branch-right")
  ts = TreeStyle()
  ts.scale = 40   
  ts.draw_guiding_lines = True
  ts.show_leaf_name = False
  ts.force_topology = True
  ts.show_scale = False
  t.render("img_faces.png", h=picture_size*0.8, tree_style=ts)


#将已经画好的物种树的图和热图合并到一张图上
def combine_fig(fig_name, picture_size):

  #先通过树图的大小计算整张图片的面积
  treepic = Image.open("img_faces.png")
  treepic_size = treepic.size
  combine_fig_size = treepic_size[0] + treepic_size[1]

  #创建画布
  combine = Image.new("RGB", (combine_fig_size, combine_fig_size), "#FFFFFF")

  #将treepic粘贴过来
  combine.paste(treepic, (int(picture_size*0.01), int(picture_size*0.01)))
  treepic_rotate = treepic.rotate(90, expand = 1) #treepic旋转90度后粘贴
  combine.paste(treepic_rotate, (treepic_size[0] - int(picture_size*0.01), 
                int(treepic_size[1]) - int(picture_size*0.01)))

  #讲hotpic粘贴过来
  hotpic = Image.open("hotmap.png")
  hotpic.thumbnail((treepic_size[1], treepic_size[1]))
  combine.paste(hotpic, (treepic_size[0] - int(picture_size*0.01), 
                 int(picture_size*0.01)))

  #保存图片
  combine.save(fig_name + ".png")

  os.remove("hotmap.png")
  os.remove("img_faces.png")


#主程序
def main():
  #解析参数
  parser = argparse.ArgumentParser(description="Options for visual_hyde.py", 
                                    add_help=True)
  parser = argparse.ArgumentParser(formatter_class=
                                   argparse.RawDescriptionHelpFormatter, 
                                    description=Description)
  required = parser.add_argument_group("Required arguments")
  required.add_argument('-i', '--infile', action="store", metavar='\b', 
                                    type=str, required=True, 
                                    help="hyde output file")  
  required.add_argument('-t', '--treefile', action="store", metavar='\b', 
                                    type=str, required=True, 
                                    help="species tree file")
  additional = parser.add_argument_group("Additional arguments")
  additional.add_argument('-n', '--node', action="store_true", default=False, 
                          help='''Node model, stack up all heatmaps for each 
                          monophyletic clade respectively, only the squares in 
                          all stacked heatmaps were light, the corresponding 
                          square will be light''')                                    
  additional.add_argument('-c', '--preclade', action="store", metavar='\b', 
                          type=str, help='''Name of predefinition clade file, 
                          different clade will have different colors on the 
                          heatmap, if not specify this file, script will 
                          defines up to five leaves as a clade automatically''')
  additional.add_argument('-l', '--leaves', action="store", metavar='\b', 
                          type=str, help='''(if not node model) Specify the name 
                          of a leaf and only draw the heatmap for it. if not 
                          specify this parameter, script will draw the heatmap
                          for all the leaves''')
  additional.add_argument('-s', '--picturesize', action="store", type=float, 
                          default=4000, metavar='\b', 
                          help="Size of the output figure, default = 4000")
  additional.add_argument('-z', '--zscore', action="store", type=float, 
                          default=3, metavar='\b', 
                          help='''threshold of Z-score of hyde output, 
                          default = 3''')


  args                           = parser.parse_args()
  csv_file_name                  = args.infile
  tree_file                      = args.treefile
  node_model                     = args.node
  Predefined_clade_file          = args.preclade
  hypothesis_hybrid_species      = args.leaves  
  picture_size                   = args.picturesize
  zscore                         = args.zscore
   
  #检查树是否置根，且外群只有一个
  input_tree = Tree(tree_file)
  if len(input_tree.children) == 3:
    print("The input tree are not rooted, script end")
    sys.exit(0)
  else:
    judge = False
    for each_node in input_tree.children:
      if each_node.is_leaf():
        judge = True
    if judge:
      pass
    else:
      print("The input tree has more than one outgroup, script end")
      sys.exit(0)      

 
  #检查树中的物种名称是否和hyde软件输出结果中的物种名称一一对应
  leaves_name_in_hyde_output = []
  with open(csv_file_name, "r") as read_file:
    for each_line in read_file.readlines()[1:]:
      if each_line.split("\t")[0] not in leaves_name_in_hyde_output:
        leaves_name_in_hyde_output.append(each_line.split("\t")[0])
      if each_line.split("\t")[1] not in leaves_name_in_hyde_output:
        leaves_name_in_hyde_output.append(each_line.split("\t")[1])
      if each_line.split("\t")[2] not in leaves_name_in_hyde_output:
        leaves_name_in_hyde_output.append(each_line.split("\t")[2])
  leaves_name_in_hyde_output.sort()
  leaves_name_in_species_tree = []
  if (len(input_tree.children[0].get_leaf_names()) > 
      len(input_tree.children[1].get_leaf_names())):
    ingroup_list = input_tree.children[0].get_leaf_names()
  else:
    ingroup_list = input_tree.children[1].get_leaf_names()
  for each_line in ingroup_list:
    if each_line not in leaves_name_in_species_tree:
      leaves_name_in_species_tree.append(each_line)
  leaves_name_in_species_tree.sort()
  if leaves_name_in_hyde_output == leaves_name_in_species_tree:
    pass
  else:
    for each_name in leaves_name_in_hyde_output:
      if each_name not in leaves_name_in_species_tree:
        print(each_name + " in hyde output but not in species tree")
    for each_name in leaves_name_in_species_tree:
      if each_name not in leaves_name_in_hyde_output:
        print(each_name + " in species tree but not in hyde output")
    print('''The sample names of ingroups in input tree and sample names of hyde file are not uniform, script end''')
    sys.exit(0)


  #检查用户是否预先设定了分支，如果没有则计算一个，这些分支会在树中具有不同的颜色
  if Predefined_clade_file:
    pass
  else:
    print('''No preclade file specified, script will defines up some clade automatically''')
    make_predefined_clade_file(tree_file)
  Predefined_clade_file = "Predefined_clade.txt"


  #开始绘制热图 
  if node_model:  #使用节点模式运行，计算各个节点的杂交情况
    print('''Run in node model, the heatmap shows the hybridization events that common to all the samples after the node''')
    #获取输入树的一些信息
    t, Label, clade_file, name_len = parse_tree(tree_file, 
                                                Predefined_clade_file)
    #将所有sample的heatmap表格存入到一个字典中，方便之后的步骤调用
    hyde_output_array_dict = make_hyde_output_array_dict(Label, csv_file_name, 
                                                         zscore)
    #开始遍历系统发育树
    node_num = 0
    for each_node in t.traverse("postorder"):
      if each_node.is_leaf():
        pass
      else:
        node_num = node_num + 1
        print("Drawing figure of node " + str(node_num) + " " + str(node_num) + "/" + str(len(t.get_leaf_names())-1))
        #将某个分支所有sample的热图叠加成一张
        hyde_output_array = find_common_hybrid_in_nodes(hyde_output_array_dict, 
                                                        each_node, Label)
        #绘制叠加好的热图
        draw_hotmap(str(node_num), hyde_output_array)
        #绘制热图旁边儿的树
        draw_tree(tree_file, node_num, each_node.get_leaf_names(), clade_file, 
                  name_len, picture_size)
        #组合热图和树
        combine_fig(str(node_num), picture_size)


  else:   #使用非节点模式，为每个sample都画一张热图
    #画图
    def run(hypothesis_hybrid_species, Predefined_clade_file, fig_num):
      t, Label, clade_file, name_len = parse_tree(tree_file, 
                                                  Predefined_clade_file)
      hyde_output_array = make_hotmap_table_gamma(Label, 
                              hypothesis_hybrid_species, csv_file_name, zscore)
      draw_hotmap(fig_num + hypothesis_hybrid_species, hyde_output_array)
      draw_tree(tree_file, "", hypothesis_hybrid_species, clade_file, 
                name_len, picture_size)
      combine_fig((fig_num + hypothesis_hybrid_species), picture_size)
    
    if hypothesis_hybrid_species: #如果用户指定了-l参数，则只画该sample
      print("Run in leaves model, the script will draw the heatmap for sample: " + hypothesis_hybrid_species)
      run(hypothesis_hybrid_species, Predefined_clade_file, "")
    else: #如果用户没有指定-l参数，则画所有sample
      print('''Run in leaves model and no leaf name specified, the script will draw the heatmap for each sample''')
      t = Tree(tree_file)
      n = 0
      for hypothesis_hybrid_species in t.get_leaf_names():
        n = n + 1
        print("Start drawing figure of sample: " + hypothesis_hybrid_species + " " + str(n) + "/" + str(len(t.get_leaf_names())))
        run(hypothesis_hybrid_species, Predefined_clade_file, str(n))        

if __name__ == "__main__":
  main()  