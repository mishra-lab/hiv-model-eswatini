export PYTHONPATH := .:$$PYTHONPATH
do.model:
	python3.8 model/main.py && echo [OK]
do.params.distr:
	python3.8 params/distr.py && echo [OK]
do.params.fsw:
	Rscript params/fsw/main.r && echo [OK]
do.params.np-gp:
	Rscript params/np-gp/main.r && echo [OK]
do.debug.py:
	python3.8 .debug/ta-save.py && echo [OK]
do.debug.r:
	Rscript .debug/duration.r && echo [OK]
do.profile: # TODO: make this its own script
	kernprof -lv -u 1 model/main.py && rm main.py.lprof && echo [OK]