PYTHON ?= python

.PHONY: verify verify-class5 verify-class6 verify-generalization verify-ml paper note note-generalization docs clean

verify: verify-class5 verify-class6 verify-generalization verify-ml

verify-class5:
	$(PYTHON) scripts/class5/verify_coherent_symbolic.py
	$(PYTHON) scripts/class5/verify_coherent_numeric.py
	$(PYTHON) scripts/class5/verify_rmatrix_comparison.py

# The existing physics CI invokes this target; include the common correction checks.
verify-class6: verify-generalization
	$(PYTHON) scripts/class6/verify_coherent_symbolic.py
	$(PYTHON) scripts/class6/verify_coherent_numeric.py

verify-generalization:
	$(PYTHON) scripts/generalization/verify_shared_site_identity.py
	$(PYTHON) scripts/generalization/verify_shared_site_proof.py

verify-ml:
	$(PYTHON) scripts/ml/verify_conventions.py
	$(PYTHON) scripts/ml/class5_ml_recovery.py
	$(PYTHON) scripts/ml/class6_ml_recovery.py

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

note:
	cd notes && latexmk -lualatex -interaction=nonstopmode -halt-on-error full_research_note_ja.tex

note-generalization:
	cd notes && latexmk -lualatex -interaction=nonstopmode -halt-on-error shared_site_generalization_ja.tex

docs: paper note

clean:
	cd paper && latexmk -C main.tex || true
	cd notes && latexmk -C full_research_note_ja.tex || true
	cd notes && latexmk -C shared_site_generalization_ja.tex || true
