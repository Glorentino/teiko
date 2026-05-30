.PHONY: setup pipeline dashboard

setup:
	pip install -r teiko_proj/requirements.txt

pipeline:
	python load_data.py
	python teiko_proj/analysis.py

dashboard:
	python teiko_proj/dashboard.py
