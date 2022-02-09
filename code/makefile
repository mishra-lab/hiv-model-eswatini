export PYTHONPATH := .:$$PYTHONPATH
do.model:
	python3.8 model/main.py && echo [OK]
do.params.distr:
	python3.8 params/distr.py && echo [OK]
do.obj.1:
	Rscript post/main_obj_1.r && echo [OK]
do.obj.2:
	Rscript post/main_obj_2.r && echo [OK]
do.params.fsw:
	Rscript params/fsw/main.r && echo [OK]
do.params.np-gp:
	Rscript params/np-gp/main.r && echo [OK]
do.sankey:
	Rscript post/sankey.r && echo [OK]