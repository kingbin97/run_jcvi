#!/bin/bash
####################################################################################################
### 使用说明                                                                                     ###
####################################################################################################
usage() {
	echo -e "\e[32mUsage: run_jcvi.sh <-q query prefix> <-t target prefix> [-j threader] [-s seq type] [-m map soft] [-d dist value] [-c cscore] [-w width] [-e height]\e[0m"
	echo -e "\e[32m       -q          字面意思                                            \e[0m"
	echo -e "\e[32m       -t          字面意思                                            \e[0m"
	echo -e "\e[32m       -j          字面意思                           [Default: 40]\e[0m"
	echo -e "\e[32m       -s          nucl | prot                        [Default: prot]\e[0m"
	echo -e "\e[32m       -m          diamond | blastp | last | blastn   [Default: diamond]\e[0m"
	echo -e "\e[32m       -x          比对的最大匹配                     [Default: 25]\e[0m"
	echo -e "\e[32m       -d          dist value                         [Default: 200]\e[0m"
	echo -e "\e[32m       -c          cscore value                       [Default: 0.9]\e[0m"
	echo -e "\e[32m       -w          dotplot width                      [Default: 15]\e[0m"
	echo -e "\e[32m       -e          dotplot height                     [Default: 15]\e[0m"
}
if [ $# -eq 0 ] || [ "$1" = "-h" ]; then
	usage
	exit 0
elif [ $# -lt 2 ]; then
	echo "Error: Missing required arguments."
	usage
	exit 1
fi
####################################################################################################
### 传参                                                                                         ###
####################################################################################################
while getopts ":q:t:j:s:m:d:c:w:e:x:" opt; do
	case $opt in
		q) A_prefix="$OPTARG";;
		t) B_prefix="$OPTARG";;
		j) THREADER="$OPTARG";;
		s) SEQ_TYPE="$OPTARG";;
		m) MAP_SOFT="$OPTARG";;
		d) DIST_VALUE="$OPTARG";;
		c) CSCORE_VALUE="$OPTARG";;
		w) WIDTH_VALUE="$OPTARG";;
		e) HEIGHT_VALUE="$OPTARG";;
		x) MAX_MAP="$OPTARG";;
		\?) echo "Invalid option -$OPTARG" >&2;;
	esac
done
####################################################################################################
# 输入文件(可以绝对路径，也可以相对路径，提取的是文件后缀之前的所有字符串)
#A=$1
#B=$2
#A_prefix=$(basename "$A" .bed)
#B_prefix=$(basename "$B" .bed)
#A_prefix=$1
#B_prefix=$2
####################################################################################################
### 线程数，空值默认40                                                                           ###
####################################################################################################
if [[ -z $THREADER ]]; then
	THREADER=40
fi
####################################################################################################
### 比对类型，空值默认 prot                                                                      ###
####################################################################################################
if [[ -z $SEQ_TYPE ]]; then
	SEQ_TYPE="prot"
fi
####################################################################################################
### 比对软件，空值默认diamond                                                                    ###
####################################################################################################
if [[ -z $MAP_SOFT ]]; then
	MAP_SOFT="diamond"
fi
####################################################################################################
### dist值，空值默认200                                                                          ###
####################################################################################################
if [[ -z $DIST_VALUE ]]; then
	DIST_VALUE=200
fi
####################################################################################################
### CSCORE 值，空值默认0.9                                                                       ###
####################################################################################################
if [[ -z $CSCORE_VALUE ]]; then
	CSCORE_VALUE=0.9
fi
####################################################################################################
### WIDTH 值，空值默认13                                                                       ###
####################################################################################################
if [[ -z $WIDTH_VALUE ]]; then
	WIDTH_VALUE=15
fi
####################################################################################################
### HEIGHT 值，空值默认 HEIGHT                                                                   ###
####################################################################################################
if [[ -z $HEIGHT_VALUE ]]; then
	HEIGHT_VALUE=$WIDTH_VALUE
fi
####################################################################################################
### 最大比对值 值，空值默认1                                                                     ###
####################################################################################################
if [[ -z $MAX_MAP ]]; then
	MAX_MAP=25
fi
####################################################################################################
### 输出参数                                                                                     ###
####################################################################################################
echo -e "\e[32m###############################################################################################################\e[0m"
echo -e "\e[32m########## 比对序列: $A_prefix                                      ##########\e[0m"
echo -e "\e[32m########## 目标序列: $B_prefix                                      ##########\e[0m"
echo -e "\e[32m########## 线程数为： $THREADER                                     ##########\e[0m"
echo -e "\e[32m########## 比对类型为: $SEQ_TYPE                                    ##########\e[0m"
echo -e "\e[32m########## 比对软件为: $MAP_SOFT                                    ##########\e[0m"
echo -e "\e[32m########## dist value: $DIST_VALUE                                  ##########\e[0m"
echo -e "\e[32m########## cscore value: $CSCORE_VALUE                              ##########\e[0m"
echo -e "\e[32m###############################################################################################################\e[0m"
####################################################################################################
### 选择比对软件模式，默认prot                                                                   ###
####################################################################################################
if [ $SEQ_TYPE == nucl ]; then
	echo "核苷酸"
	echo "Blastn比对模式"
elif [ $SEQ_TYPE == "prot" ]; then
	if [ $MAP_SOFT == "blastp" ]; then
		echo "Blastp比对模式"
		echo "      "
	elif [ $MAP_SOFT == "last" ]; then
		echo "last模式"
		echo "      "
	elif [ $MAP_SOFT == "diamond" ]; then
		echo "diamond比对模式"
		echo "      "
	else
		echo "是prot，请检查比对模式"
		echo "      "
		exit 1
	fi
else
	echo "请检查数据"
	exit 1
fi
####################################################################################################
### function                                                                                     ###
####################################################################################################
function get_cds_A() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m get ${A_prefix} cds\e[0m"
	gffread -g ${A_prefix}.fasta -w ${A_prefix}.raw.cds ${A_prefix}.gff3 #####提取raw cds
	bioawk -c fastx 'length($seq) > 90{ print ">"$name; print $seq }' ${A_prefix}.raw.cds > ${A_prefix}.raw.90bp.cds #####挑选大于90bp
	removeRedundantProteins_2.py -i ${A_prefix}.raw.90bp.cds -o ${A_prefix}.cds #####去冗余，挑选最长转录本
	#rm ${A_prefix}.raw.cds ${A_prefix}.raw.90bp.cds #####删除中间文件 
}
function get_cds_B() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m get ${B_prefix} cds\e[0m"
	gffread -g ${B_prefix}.fasta -w ${B_prefix}.raw.cds ${B_prefix}.gff3
	bioawk -c fastx 'length($seq) > 90{ print ">"$name; print $seq }' ${B_prefix}.raw.cds > ${B_prefix}.raw.90bp.cds
	removeRedundantProteins_2.py -i ${B_prefix}.raw.90bp.cds -o ${B_prefix}.cds
	#rm ${B_prefix}.raw.cds ${B_prefix}.raw.90bp.cds
}
function get_pep_A() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m get${A_prefix} pep\e[0m"
	gffread -g ${A_prefix}.fasta -V -y ${A_prefix}.raw.pep ${A_prefix}.gff3
	bioawk -c fastx 'length($seq) > 30{ print ">"$name; print $seq }' ${A_prefix}.raw.pep > ${A_prefix}.raw.30aa.pep
	removeRedundantProteins_2.py -i ${A_prefix}.raw.30aa.pep -o ${A_prefix}.pep
	#rm ${A_prefix}.raw.pep ${A_prefix}.raw.30aa.pep
}
function get_pep_B() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m get${B_prefix} pep\e[0m"
	gffread -g ${B_prefix}.fasta -V -y ${B_prefix}.raw.pep ${B_prefix}.gff3
	bioawk -c fastx 'length($seq) > 30{ print ">"$name; print $seq }' ${B_prefix}.raw.pep > ${B_prefix}.raw.30aa.pep
	removeRedundantProteins_2.py -i ${B_prefix}.raw.30aa.pep -o ${B_prefix}.pep
	#rm ${B_prefix}.raw.pep ${B_prefix}.raw.30aa.pep
}
function get_bed_A() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m get${A_prefix} bed\e[0m"
	python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only ${A_prefix}.gff3 -o ${A_prefix}_primary_only.bed
	python -m jcvi.formats.bed uniq ${A_prefix}_primary_only.bed
	cp ${A_prefix}_primary_only.uniq.bed ${A_prefix}.bed
}
function get_bed_B() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m get${B_prefix} bed\e[0m"
	#python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only ${B_prefix}.gff3 -o ${B_prefix}.bed
	python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only ${B_prefix}.gff3 -o ${B_prefix}_primary_only.bed
	python -m jcvi.formats.bed uniq ${B_prefix}_primary_only.bed
	cp ${B_prefix}_primary_only.uniq.bed ${B_prefix}.bed
}
function pairwise_synteny_01() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m##### pairwise_synteny_01\e[0m"
	###01.共线性对可视化
	python -m jcvi.compara.catalog ortholog ${A_prefix} ${B_prefix} --no_strip_names --nochpf --nostdpf --skipempty --nostdp --dist=${DIST_VALUE} --cscore=${CSCORE_VALUE} --dbtype=${SEQ_TYPE} --cpus=${THREADER} --self_remove=${self_remove} #--no_dotplot
	### dotplot
	python -m jcvi.graphics.dotplot ${A_prefix}.${B_prefix}.lifted.anchors --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --font=Arial --format=pdf -o ${A_prefix}.${B_prefix}.lifted.anchors.dotplot.pdf --minfont=3 --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE} --nochpf --nostdpf --dpi=600 --nosort --skipempty #--notex 
	python -m jcvi.graphics.dotplot ${A_prefix}.${B_prefix}.lifted.anchors --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --font=Arial --format=png -o ${A_prefix}.${B_prefix}.lifted.anchors.dotplot.png --minfont=3 --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE} --nochpf --nostdpf --dpi=600 --nosort --skipempty #--notex 
	python -m jcvi.graphics.dotplot ${A_prefix}.${B_prefix}.lifted.anchors --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --font=Arial --format=svg -o ${A_prefix}.${B_prefix}.lifted.anchors.dotplot.svg --minfont=3 --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE} --nochpf --nostdpf --dpi=600 --nosort --skipempty #--notex 
	### 彩图
	python -m jcvi.graphics.dotplot ${A_prefix}.${B_prefix}.lifted.anchors --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --font=Arial --format=pdf -o ${A_prefix}.${B_prefix}.lifted.anchors.color.pdf --minfont=3 --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE} --nochpf --nostdpf --theme 1 --style white --colororientation --dpi=600 --nosort --skipempty #--notex 
	python -m jcvi.graphics.dotplot ${A_prefix}.${B_prefix}.lifted.anchors --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --font=Arial --format=png -o ${A_prefix}.${B_prefix}.lifted.anchors.color.png --minfont=3 --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE} --nochpf --nostdpf --theme 1 --style white --colororientation --dpi=600 --nosort --skipempty #--notex 
	python -m jcvi.graphics.dotplot ${A_prefix}.${B_prefix}.lifted.anchors --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --font=Arial --format=svg -o ${A_prefix}.${B_prefix}.lifted.anchors.color.svg --minfont=3 --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE} --nochpf --nostdpf --theme 1 --style white --colororientation --dpi=600 --nosort --skipempty #--notex 
	### 染色体坐标图
	python -m jcvi.graphics.blastplot ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.filtered --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE}  --dpi=600 --format=png
	mv ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.png ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.filtered.chrome.png
	python -m jcvi.graphics.blastplot ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.filtered --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE}  --dpi=600 --format=pdf
	mv ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.pdf ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.filtered.chrome.pdf
	python -m jcvi.graphics.blastplot ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.filtered --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE}  --dpi=600 --format=svg
	mv ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.svg ${A_prefix}.${B_prefix}.last.P${self_remove}L0.inverse.filtered.chrome.svg
}
function pairwise_synteny_RBH_02() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m##### pairwise_synteny_RBH_02\e[0m"
	####02.pairwise synteny RBH
	python -m jcvi.compara.catalog ortholog ${A_prefix} ${B_prefix} --full --no_strip_names --nochpf --nostdpf --skipempty --nostdp --dist=${DIST_VALUE} --cscore=${CSCORE_VALUE} --dbtype=${SEQ_TYPE} --cpus=${THREADER} --self_remove=${self_remove} --no_dotplot
}
function pairwise_synteny_depth_03() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m##### pairwise_synteny_depth_03\e[0m"
	####03.synteny depth
	python -m jcvi.compara.synteny depth ${A_prefix}.${B_prefix}.anchors --histogram --depthfile=${A_prefix}.${B_prefix}.anchors.depth.txt
	#python -m jcvi.compara.synteny depth ${A_prefix}.${B_prefix}.anchors --histogram
}
function blastfilter_04() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m##### blastfilter_04\e[0m"
	####04.blast结果筛选，如果需要，可以把filter的结果作为共线性分析输入的初始文件
	python -m jcvi.compara.blastfilter ${A_prefix}.${B_prefix}.last --qbed ${A_prefix}.bed --sbed ${B_prefix}.bed --cscore=${CSCORE_VALUE} --no_strip_names
	###染色体坐标图
	#python -m jcvi.graphics.blastplot ${A_prefix}.${B_prefix}.last.filtered --sbed ${B_prefix}.bed --qbed ${A_prefix}.bed --figsize=${WIDTH_VALUE}x${HEIGHT_VALUE} --dpi=600 --font=Arial --style=white --notex
	#mv ${A_prefix}.${B_prefix}.last.pdf ${A_prefix}.${B_prefix}.last.filtered.chrome.pdf
}
function bezier_curve_graph() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo -e "\e[31m##### bezier_curve_graph\e[0m"
	## 01.首先生成.simpleple文件
	python -m jcvi.compara.synteny screen --minspan=8 --simple ${A_prefix}.${B_prefix}.anchors ${A_prefix}.${B_prefix}.anchors.new
	## 02.生成设置需要展示的染色体号的seqids文件，
	cut -f1 ${A_prefix}.bed |sort |uniq |sed ':a;N;$!ba;s/\n/,/g' > ${A_prefix}.seqid.txt
	cut -f1 ${B_prefix}.bed |sort |uniq |sed ':a;N;$!ba;s/\n/,/g' >${B_prefix}.seqid.txt
	A_line_count=$(cut -f1 ${A_prefix}.bed | sort | uniq | wc -l) # 获取A文件行数
	B_line_count=$(cut -f1 ${B_prefix}.bed | sort | uniq | wc -l) # 获取B文件行数
	max_value=$((A_line_count > B_line_count ? A_line_count : B_line_count)) # 获取最大值
	echo "$max_value"
	if (( max_value >= 1 && max_value <= 30 )); then
		bezier_width=8
	elif (( max_value >= 31 && max_value <= 60 )); then
		bezier_width=12
	else
		bezier_width=16
	fi
	cat ${A_prefix}.seqid.txt ${B_prefix}.seqid.txt >${A_prefix}.${B_prefix}.seqid.txt
	#rm  ${A_prefix}.seqid.txt ${B_prefix}.seqid.txt
	## 03.生成设置颜色，长宽等的layout文件，
	echo -e "# y, xstart, xend, rotation, color, label, va, bed\n 0.7, .25, .75,0, , ${A_prefix}, top, ./${A_prefix}.bed\n 0.3, .25, .75,0, , ${B_prefix}, bottom, ./${B_prefix}.bed\n# edges\ne, 0, 1, ./${A_prefix}.${B_prefix}.anchors.simple" > ${A_prefix}.${B_prefix}.layout.txt
	## 04.运行代码
	python -m jcvi.graphics.karyotype ${A_prefix}.${B_prefix}.seqid.txt ${A_prefix}.${B_prefix}.layout.txt --font=Arial --notex --format=svg --dpi=600 --figsize=${bezier_width}x${bezier_width}
}
function blastn_map() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo "blastn_map"
	####建索引
	makeblastdb -in ${B_prefix}.cds -dbtype nucl -title ${B_prefix} -out ${B_prefix}
	####比对 #-max_hsps 1 -num_alignments 20
	blastn -db ${B_prefix} -query ${A_prefix}.cds -out ${A_prefix}.${B_prefix}.last -evalue 1e-5 -outfmt 6 -num_threads ${THREADER} -max_hsps 1 -num_alignments ${MAX_MAP}
}
function blastp_map() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo "blastp_map"
	####建索引
	diamond makedb --in ${B_prefix}.pep --db ${B_prefix}.pep
	####比对 #-max_hsps 1 -num_alignments 20
	blastp -db ${B_prefix} -query ${A_prefix}.pep -out ${A_prefix}.${B_prefix}.last -evalue 1e-5 -outfmt 6 -num_threads ${THREADER} -max_hsps 1 -num_alignments ${MAX_MAP}
}
function diamond_map() {
	echo -e "\e[32m###############################################################################################################\e[0m"
	echo "diamond_map"
	####建索引
	diamond makedb --in ${B_prefix}.pep --db ${B_prefix}.pep
	####比对 --max-hsps 1 --max-target-seqs 25 --no-self-hits
	diamond blastp --query ${A_prefix}.pep --db ${B_prefix}.pep --evalue 1e-5 --very-sensitive --outfmt 6 --threads ${THREADER} --quiet --header --out ${A_prefix}.${B_prefix}.last --max-hsps 1 --max-target-seqs ${MAX_MAP} --no-self-hits 
	rm ${B_prefix}.pep.dmnd
}
####################################################################################################
### 检测是否有bed文件和last文件，如果有则直接执行共线性分析                                      ###
####################################################################################################

####################################################################################################
### 检测是否有A.bed文件                                                                          ###
####################################################################################################
if [ ! -e ${A_prefix}.bed ] ; then  ###没有${A_prefix}.bed文件
	if [ -e ${A_prefix}.gff3 ]; then ###有${A_prefix}.gff文件
		get_bed_A
		echo "get_bed_A"
	else
		echo "请检查${A_prefix}的bed文件，或者gff3文件"
		echo "      "
		exit 1
	fi
else
	echo "有${A_prefix}.bed文件"
fi
####################################################################################################
### 检测是否有B.bed文件                                                                          ###
####################################################################################################
if [ ! -e ${B_prefix}.bed ] ; then
	if [ -e ${B_prefix}.gff3 ]; then ###有${B_prefix}.gff文件
		get_bed_B
		echo "get_bed_B"
		echo "      "
	else
		echo "请检查${B_prefix}的bed文件，或者gff3文件"
		echo "      "
		exit 1
	fi
else
	echo "有${B_prefix}.bed文件"
fi
####################################################################################################
### 检测是否有A.B.last文件                                                                       ###
####################################################################################################
if [ -e ${A_prefix}.${B_prefix}.last ]; then
	echo "检测到 ${A_prefix}.${B_prefix}.last 文件"
else
	case $SEQ_TYPE in
		prot) ##蛋白质模式
			echo "蛋白质模式"
			if [[ -e "${A_prefix}.pep" && -e "${B_prefix}.pep" ]]; then
				echo "检测到 ${A_prefix}.pep 与 ${B_prefix}.pep ,启动蛋白质比对模式"
			else
				##检测A.pep文件
				if [[ -e "${A_prefix}.pep" ]]; then
					echo "检测到${A_prefix}.pep"
				else
					echo "未检测到${A_prefix}.pep，将检测 ${A_prefix}.fasta与${A_prefix}.gff3"
					if [[ -e "${A_prefix}.fasta" && -e "${A_prefix}.gff3" ]]; then
						echo "检测到 ${A_prefix}.fasta 与 ${A_prefix}.gff3"
						get_pep_A
						echo "get {A_prefix} pep"
					else
						echo "请检查 ${A_prefix}.fasta 与 ${A_prefix}.gff3 数据"
						exit 1
					fi
				fi
				##检测B.pep文件
				if [[ -e "${B_prefix}.pep" ]]; then
					echo "      "
					echo "检测到${B_prefix}.pep"
				else
					echo "      "
					echo "未检测到${B_prefix}.pep，将检测 ${B_prefix}.fasta与${B_prefix}.gff3"
					if [[ -e "${B_prefix}.fasta" && -e "${B_prefix}.gff3" ]]; then
						echo "检测到 ${B_prefix}.fasta 与 ${B_prefix}.gff3"
						get_pep_B
						echo "get {B_prefix} pep"
						
					else
						echo "请检查 ${B_prefix}.fasta 与 ${B_prefix}.gff3 数据"
						exit 1
					fi
				fi
			fi
			##得到A.pep和B.pep后，进行比对
			if [[ $MAP_SOFT == "diamond" ]]; then
				echo "      "
				echo "diamond模式"
				diamond_map
			elif [[ $MAP_SOFT == "blastp" ]]; then
				echo "      "
				echo "blastp模式"
				blastp_map
			elif [[ $MAP_SOFT == "last" ]]; then
				echo "      "
				echo "last模式"
			else
				echo "      "
				echo "请检查比对模式"
				exit 1
			fi
		;;
		nucl)
			echo "核苷酸比对模式"
			if [[ -e "${A_prefix}.cds" && -e "${B_prefix}.cds" ]]; then
				echo "检测到 ${A_prefix}.cds 与 ${B_prefix}.cds ,启动核苷酸比对模式"
			else
				##检测A.cds
				if [[ -e "${A_prefix}.cds" ]]; then
					echo "检测到 ${A_prefix}.cds"
				else
					echo "未检测到 ${A_prefix}.cds，将检测 ${A_prefix}.fasta与${A_prefix}.gff3"
					if [[ -e "${A_prefix}.fasta" && -e "${A_prefix}.gff3" ]]; then
						echo "检测到 ${A_prefix}.fasta 与 ${A_prefix}.gff3，将提取 ${A_prefix}.cds"
						get_cds_A
						echo "get {A_prefix} cds"
					else
						echo "未检测到 ${A_prefix}.fasta 与 ${A_prefix}.gff3"
						exit 1
					fi
				fi
				##检测B.cds
				if [[ -e "${B_prefix}.cds" ]]; then
					echo "检测到 ${B_prefix}.cds"
				else
					echo "未检测到 ${B_prefix}.cds，将检测 ${B_prefix}.fasta与${B_prefix}.gff3"
					if [[ -e "${B_prefix}.fasta" && -e "${B_prefix}.gff3" ]]; then
						echo "检测到 ${B_prefix}.fasta 与 ${B_prefix}.gff3，将提取 ${B_prefix}.cds"
						get_cds_B
						echo "get {B_prefix} cds"
					else
						echo "未检测到 ${B_prefix}.fasta 与 ${B_prefix}.gff3"
						exit 1
					fi
				fi
			fi
			##若上述cds准备ok，进行下一步比对
			echo "检测到 ${A_prefix}.cds 与 ${B_prefix}.cds ,启动核苷酸比对模式"
			blastn_map
			echo "blastn_map"
		;;
		*)
			echo -e "\e[32m###############################################################################################################\e[0m"
			echo "请检测输入文件或参数"
			exit 1
		;;
	esac
fi
####################################################################################################
### colinearity                                                                                  ###
####################################################################################################
if [[ -e "${A_prefix}.bed" && -e "${B_prefix}.bed" && -e "${A_prefix}.${B_prefix}.last" ]]; then
	echo "      "
	echo "检测到 bed文件 和 last文件"
	echo "      "
	self_remove=97
	pairwise_synteny_01
	echo "      "
	#if [[ -e "${A_prefix}.cds" && -e "${B_prefix}.cds" && -e "${A_prefix}.${B_prefix}.last" ]]; then
	#	echo "RBH模式"
	#	pairwise_synteny_RBH_02
	#else
	#	echo "未检测到 ${A_prefix}.cds 与 ${B_prefix}.cds , 将跳过RBH模式"
	#fi
	echo "      "
	pairwise_synteny_depth_03
	echo "      "
	#blastfilter_04
	echo "      "
	bezier_curve_graph
	mv karyotype.pdf ${A_prefix}.${B_prefix}.karyotype.pdf
	exit 0
else
	echo "      "
	echo "请检查文件"
	exit 1
fi
