cat config | while read id
do
	arr=($id)
	name1=${arr[0]}
    name2=${arr[1]}
    dir=${name1}_${name2}
    mkdir ${dir}
    cp ${name1}.uniq.bed ${dir}/${name1}.bed
	cp ${name2}.uniq.bed ${dir}/${name2}.bed
	cp ${name1}.cds ${dir}/
	cp ${name2}.cds ${dir}/
	cd ${dir}
	python -m jcvi.compara.catalog ortholog --no_strip_names ${name1} ${name2} --cscore=.99
	python -m jcvi.compara.synteny screen --minspan=30 --simple ${name1}.${name2}.anchors ${name1}.${name2}.anchors.new
	tmpstr1=`less -S ${name1}.bed | awk '{print $1}' | uniq | tr '\n' ','`
	str1=${tmpstr1%?}
	echo $str1 > seqids
	tmpstr2=`less -S ${name2}.bed | awk '{print $1}' | uniq | tr '\n' ',' `
	str2=${tmpstr2%?}
	echo $str2 >> seqids
content="
# y, xstart, xend, rotation, color, label, va,  bed\n .6,     .1,    .8,       0,      , ${name1#*_}, top, ${name1}.bed\n .4,     .1,    .8,       0,      , ${name2#*_}, top, ${name2}.bed\n# edges\ne, 0, 1, ${name1}.${name2}.anchors.simple
"

echo -e $content > layout

python -m jcvi.graphics.karyotype seqids layout
cd ../

done
