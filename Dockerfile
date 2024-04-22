
FROM python:3.9

ARG ssh_prv_key
ARG ssh_pub_key

RUN apt-get update && \
    apt-get install -y \
        git \
        openssh-server 
        # \libmysqlclient-dev

RUN echo "ls"

RUN mkdir -p /root/.ssh && \
    chmod 0700 /root/.ssh

COPY ./requirements.txt /code/requirements.txt
# See: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/githubs-ssh-key-fingerprints
COPY ./known_hosts  /root/.ssh/known_hosts

# Add the keys and set permissions
RUN echo "$ssh_prv_key" > /root/.ssh/id_rsa && \
    echo "$ssh_pub_key" > /root/.ssh/id_rsa.pub && \
    chmod 600 /root/.ssh/id_rsa && \
    chmod 600 /root/.ssh/id_rsa.pub

# Avoid cache purge by adding requirements first
ADD ./requirements.txt /code/requirements.txt
# COPY requirements.txt requirements.txt


WORKDIR /code
ENV FLASK_APP app.py
ENV FLASK_RUN_HOST 0.0.0.0
# RUN apk add --no-cache gcc musl-dev linux-headers
RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt 

# Remove SSH keys
RUN rm -rf /root/.ssh/

# Add the rest of the files
# ADD . .

COPY / /code

CMD ["flask", "run"]
