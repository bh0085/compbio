import re 

primers3 =re.findall(re.compile('[ATGUC]+'), re.sub(re.compile('[\s\-]'),'',u'''(3 ext.1: TTCTAATAC- GACTCACTA TAGGGAGGACGATGCGG), (3 ext.2: 5TCGG G CGAGTCGTCTG 3), and (3 ext. 3: 5 CTGGT CATGCGCG- GCATT TAATTCTCGGGCGAGTCGTCTG 3).'''))

primers5 = re.findall(re.compile('[AUTGC]+'), re.sub(re.compile('[\s\-]'),'',u'''(5 ext.1: 5 TCG- GGCGAGTCGTCTG 3), (5 ext.2: 5 GGGGAGAATA TTGAAATATAAATGGGAGGACGATGCGGACCG 3) and (5ext.3: 5 TTCTAAT ACGACTCACTATAGGGGAGAATATT- GAAATATAAAT 3).'''))

A9, =  re.findall(re.compile('[ATGUC]{3,}'), re.sub(re.compile('[\s\-]'),'',u'''PLA. 5GGGAGGACGAUGCGGACCGAAAAAGACCUGA CUUCUAUACUAAGUCUACGUUCCCAGACGACUCGC CCGA 3 (18).'''))

A9_dna =  re.findall(re.compile('[ATGUC]{3,}'), re.sub(re.compile('[\s\-]'),'', '''The DNA template for transcription of this aptamer had the sequence 5 TTCTAATACGACTCACTAT- AGGGAGGACGATGCG GACCGAAAAAGACCTGACT- TCTATACTAAGTCTACGTTCCCAGACGACTCGC CCGA 3.'''))

A9_3ex, A9_5ex=  re.findall(re.compile('[ATGUC]{3,}'), re.sub(re.compile('[\s\-]'),'',u'''3 extended aptamer A9
5GGGAGGACGAUGCGGACCGAAAAAGACCUGA- CUUCUAUACUAAGUCUACGUUCCCAGACGACU CGCCCGAGAAUUAAAUGCCCGCCAUGACCAG 3
2.4. Real-Time PCR
5 extended aptamer A9
5GGGAGAAUAUUGAAAUAUAAAUGGGAGGAC- GAUGCGGACCGAAAAAGACCUGACUUCUAUAC- UAAGUCUACGUUCCCAGACGACUCGCCCGA 3.'''))

print A9
print A9_3ex
print A9[::-1]
print A9_5ex[::-1]
