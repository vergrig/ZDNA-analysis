# ZDNA-analysis
Analyzing H3K9me3_ZDNA_human (commentary in Russian) 

## Анализ пиков гистоновой метки

### Скачиваем данные экспериментов
Скачиваем данные, оставляем первые 5 столбцов.

`zcat ENCFF921OTR.bed.gz | cut -f1-5 > H3K9me9_SJSA1.ENCFF921OTR.hg38.bed`

`zcat ENCFF157SWY.bed.gz | cut -f1-5 > H3K9me9_SJSA1.ENCFF157SWY.hg38.bed`

**Гистограммы длин:**
Строим их с помощью скрипта `len_hist.R`. Пиков 37100 и 42012 соответственно. 

<p float="left">
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/b9685b96d9763b3c854b21cd264da8883aa1dadd/img/len_hist.H3K9me9_SJSA1.ENCFF157SWY.hg38.png" width="450" />
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/b9685b96d9763b3c854b21cd264da8883aa1dadd/img/len_hist.H3K9me9_SJSA1.ENCFF921OTR.hg38.png" width="450" /> 
</p>

**Переводим координаты из hg38 в hg19**
Скачиваем необходимые данные о переводе, запускаем liftOver. Пиков 36888 и 41770 соответственно. 

`wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz`

`liftOver H3K9me9_SJSA1.ENCFF921OTR.hg38.bed hg38ToHg19.over.chain.gz H3K9me9_SJSA1.ENCFF921OTR.hg19.bed H3K9me9_SJSA1.ENCFF921OTR.unmapped.bed`

`liftOver H3K9me9_SJSA1.ENCFF157SWY.hg38.bed hg38ToHg19.over.chain.gz H3K9me9_SJSA1.ENCFF157SWY.hg19.bed H3K9me9_SJSA1.ENCFF157SWY.unmapped.bed`

**Гистограммы длин:** 

<p float="left">
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/len_hist.H3K9me9_SJSA1.ENCFF157SWY.hg19.png" width="450" />
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/len_hist.H3K9me9_SJSA1.ENCFF921OTR.hg19.png" width="450" /> 
</p>

### Отбросим outliers

Мы выбрали отбросить пики длиннее, чем 700, так график получается достаточно равномерным, скрипт `filter_peaks.R`. Пиков стало 31989 и 36813 соответственно.

<p float="left">
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/filter_peaks.H3K9me9_SJSA1.ENCFF157SWY.hg19.filtered.hist.png" width="450" />
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/filter_peaks.H3K9me9_SJSA1.ENCFF921OTR.hg19.filtered.hist.png" width="450" /> 
</p>

### Анализ расположения меток
Смотрим, где располагаются пики гистоновой метки относительно аннотированных генов, скрипт `chip_seeker.R`

<p float="left">
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/chip_seeker.H3K9me9_SJSA1.ENCFF157SWY.hg19.filtered.plotAnnoPie.png" width="450" /> 
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/chip_seeker.H3K9me9_SJSA1.ENCFF921OTR.hg19.filtered.plotAnnoPie.png" width="450" />  

</p>

### Объединение наборов
Объединяем два набора отфильтрованных ChIP-seq пиков с помощью утилиты bedtools merge.

`cat  *.filtered.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  H3K9me9_SJSA1.merge.hg19.bed`

![](https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/merge.PNG)

### Вторичная структура DeepZ
Скачиваем файл со вторичной стр-рой ДНК, строим распределение длин участков вторичной стр-ры ДНК, смотрим, где располагаются участки стр-ры ДНК относительно аннотированных генов.

`wget https://raw.githubusercontent.com/vanya-antonov/hse21_H3K4me3_ZDNA_human/main/data/DeepZ.bed`

<p float="left">
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/len_hist.DeepZ.png" width="450" />
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/chip_seeker.DeepZ.plotAnnoPie.png" width="450" /> 
</p>

## Анализ пересечений гистоновой метки и стр-ры ДНК

## Пересечения гистоновой меткой и стр-рами ДНК
Находим пересечения, строим гистограмму. В пересечении 881 пиков.

`bedtools intersect -a DeepZ.bed -b H3K9me9_SJSA1.merge.hg19.bed > H3K9me9_SJSA1.intersect_with_DeepZ.bed`
![]()
<p float="left">
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/len_hist.H3K9me9_SJSA1.intersect_with_DeepZ.png" width="450"/>
  <img src="https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/chip_seeker.H3K9me9_SJSA1.intersect_with_DeepZ.plotAnnoPie.png" width="450" /> 
</p>


## Визуализируем в геномном браузере
Визуализируем в геномном браузере исходные участки стр-ры ДНК, а также их пересечения с гистоновой меткой, сессия: http://genome.ucsc.edu/s/pvgrigorev/hse21_H3K9me3_ZDNA_human

![](https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/intersect1.PNG)
`chr10:103990644-103990659`

![](https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/intersect2.PNG)
`chr10:104210482-104210530`


## GO-анализ для полученных уникальных генов
Ассоциируем полученные пересечения с ближайшими генами, скрипт `ChIPpeakAnno.R`

![](https://github.com/vergrig/ZDNA-analysis/blob/d33a65e5ab5834ed5c394a1262152daf83768f68/img/pantherdb.PNG)
