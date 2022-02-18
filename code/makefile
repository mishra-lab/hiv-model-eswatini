export PYTHONPATH := .:$$PYTHONPATH
do.model:
	python3.8 model/main.py 0 && echo [OK]
do.params.distr:
	python3.8 params/distr.py && echo [OK]
do.profile:
	./profile/run.sh && echo [OK]
do.infections:
	Rscript post/infections.r && echo [OK]
do.obj.1:
	Rscript post/main_obj_1.r && echo [OK]
do.obj.2:
	Rscript post/main_obj_2.r && echo [OK]
do.params.fsw:
	Rscript params/fsw/main.r && echo [OK]
do.params.np-gp:
	Rscript params/np-gp/main.r && echo [OK]
do.model.batch:
	python3.8 model/main.py 0 && \
	python3.8 model/main.py 1 && \
	python3.8 model/main.py 2 && \
	python3.8 model/main.py 3 && \
	python3.8 model/main.py 4 && \
	python3.8 model/main.py 5 && \
	python3.8 model/main.py 6 && \
	python3.8 model/main.py 7 && \
	python3.8 model/main.py 8 && \
	python3.8 model/main.py 9 && \
	echo [OK]
