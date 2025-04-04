FROM python:3.12
WORKDIR /project
COPY pyproject.toml /project/

RUN pip install .
COPY dagster_airglow/ /project/dagster_airglow/
RUN pip install -e .
