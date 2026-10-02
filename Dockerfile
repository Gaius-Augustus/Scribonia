# Scribonia: genetic-code stop-set classifier for genomes and bins.
#
#   docker build -t gaiusaugustus/scribonia:0.1.1 .      (tag = scribonia.__version__)
#   singularity build scribonia.sif docker://gaiusaugustus/scribonia:0.1.1
#
# Feature extraction and the rule-based fallback need numpy only; the shipped
# model needs scikit-learn and joblib, which also make `scribonia train` work.
FROM python:3.12-slim

LABEL org.opencontainers.image.source="https://github.com/Gaius-Augustus/Scribonia" \
      org.opencontainers.image.version="0.1.1"

WORKDIR /opt/scribonia
COPY pyproject.toml README.md LICENSE ./
COPY scribonia ./scribonia
# The shipped model is a scikit-learn pickle: keep the version it was trained
# with (scribonia_v3).
RUN pip install --no-cache-dir ".[train]" "scikit-learn==1.9.1"
# The build fails if the installed package cannot load its model, or loads
# it with another scikit-learn version than it was trained with.
RUN cd / && python -c "import warnings; \
from sklearn.exceptions import InconsistentVersionWarning; \
warnings.simplefilter('error', InconsistentVersionWarning); \
from scribonia.model import load_model; assert load_model() is not None"

ENTRYPOINT ["scribonia"]
CMD ["--help"]
