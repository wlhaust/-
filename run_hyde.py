'''
运行hyde
输入串联起来的fasta序列
fasta序列可以是任何名称
'''
import argparse
import sys
import os


#fasta序列转化为phy文件，并且返回fasta的物种数目以及位点数目
def fasta2phy(fasta_file):
    with open(fasta_file, "r") as read_file:
        species_num = 0
        species_len = 0
        for each_line in read_file:
            if len(each_line) > 2:
                if each_line[0] == ">":
                    species_num = species_num + 1
                else:
                    species_len = len(each_line.replace("\n", ""))
    with open(fasta_file.replace(".fasta", "") + ".phy", "a") as write_file:
        with open(fasta_file, "r") as read_file:
            write_file.write(" " + str(species_num) + " " + str(species_len) + "\n")
            for each_line in read_file:
                if len(each_line) > 2:
                    if each_line[0] == ">":
                        write_file.write(each_line.replace("\n","").replace(">","") + " ")
                    else:
                        write_file.write(each_line.replace("\n","") + "\n")
    return str(species_num), str(species_len)
       
#写出map.txt文件
def make_mapfile(fasta_file, outgroup):
    with open("map.txt", "a") as write_file:
        with open(fasta_file, "r") as read_file:
            for each_line in read_file:
                if len(each_line) > 2:
                    if each_line[0] == ">":
                        if each_line.replace(">", "").replace("\n", "") == outgroup:
                            write_file.write(each_line.replace(">", "").replace("\n", "") + " " + "out\n")
                        else:
                            write_file.write(each_line.replace(">", "").replace("\n", "") + " " + each_line.replace(">", "").replace("\n", "") + "\n")

'''
if __name__ == "__main__":
    """
    Runs the script.
    """
    fasta_file_list = get_file_list()
    for each_fasta in fasta_file_list:
        species_num, species_len = fasta2phy(each_fasta)
        make_mapfile(each_fasta)
        infile = each_fasta.replace(".fasta", "") + ".phy"
        mapfile = "map.txt"
        outgroup = "out"
        nind = species_num
        ntaxa = species_num
        nsites = species_len
        os.system("python hyde.py -i " + infile + " -j 6" + " --ignore_amb_sites -m " + mapfile + " -o " + outgroup + " -n " + nind + " -t " + ntaxa + " -s " + nsites + " --prefix " + each_fasta.replace(".fasta", ""))
'''

def main():
    #解析参数
    #参数分为必须参数(required)和可选参数(additional)
    parser = argparse.ArgumentParser(description="Options for visual_hyde.py", add_help=True)
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-i', '--infile', action="store", metavar='\b', type=str, required=True, help="Name of the hyde input fasta file")  
    required.add_argument('-o', '--outgroup', action="store", metavar='\b', type=str, required=True, help="Name of the outgroup")


    args                           = parser.parse_args()
    input_fasta                    = args.infile
    outgroup                       = args.outgroup
    
    species_num, species_len = fasta2phy(input_fasta)
    make_mapfile(input_fasta, outgroup)
    infile = input_fasta.replace(".fasta", "") + ".phy"
    mapfile = "map.txt"
    outgroup = "out"
    nind = species_num
    ntaxa = species_num
    nsites = species_len
    os.system("python hyde.py -i " + infile + " -j 60" + " --ignore_amb_sites -m " + mapfile + " -o " + outgroup + " -n " + nind + " -t " + ntaxa + " -s " + nsites + " --prefix " + input_fasta.replace(".fasta", ""))


main()
os.remove("map.txt")