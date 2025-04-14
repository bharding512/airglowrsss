FROM python:3.12
WORKDIR /project
COPY pyproject.toml /project/

RUN pip install .
COPY src/airglow/ /project/airglow/
RUN pip install -e .
