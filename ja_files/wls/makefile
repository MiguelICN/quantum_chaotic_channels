.PHONY: run_wls

run_wls:
	@BASENAME=$$(basename $(file) .wls); \
	TIMESTAMP=$$(date +"%y%m%d_%H%M"); \
	LOGFILE=log/$${BASENAME}_$${TIMESTAMP}; \
	mkdir -p log; \
	(nohup ./$(file) > $${LOGFILE}.log 2>&1 & echo $$! > $${LOGFILE}.pid); \
	echo "Log output guardado en $${LOGFILE}.log \n\
	PID guardado en $${LOGFILE}.pid"
