#!/bin/bash

for files in $(ls *.eps --color=none);
do
	caption=$(echo $files | sed 's/_/\\\_/g');
	#echo $caption;
	echo '\begin{figure}';
	echo $'\t''\epspic{0.8}{'$files'}';
	echo $'\t''\caption['${caption%.eps}'\@.]{\label{fig:'${files%.eps}'}'${caption%.eps}'\@.}';
	echo '\end{figure}';
	echo '\clearpage';
	echo '';
done
