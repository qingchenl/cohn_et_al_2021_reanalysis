proc import datafile="/home/qling0/DataFolder/table_S3.xlsx" out=work.cohn
	DBMS=XLSX;
	getnames = yes;
	sheet = "A";
run;

proc print data=cohn;
run;

proc contents data=cohn;
run;


proc transpose data=cohn out=infections;
	by week;
	
run;

data infections_clean;
	keep week vaccine infected count;
	set infections;
	vaccine = prxchange('s/ (?:\+|\-)//', -1, _LABEL_);
	if find(_LABEL_, '+') then
		infected = 1;
	else
		infected = 0;
	count = col1;
	put _LABEL_ infected;
run;

proc sort data=infections_clean;
	by vaccine;
run;

proc print data = infections_clean;
run;

/*
ods graphics on;
ods select survivalplot(persist);
*/
proc lifetest data=infections_clean method=KM;
	/*by vaccine;*/
	time week*infected(0);
	strata vaccine;
	freq count;
run;
/*ods graphics off;*/

