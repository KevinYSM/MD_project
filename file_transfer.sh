#!/bin/bash
sudo scp -i /home/ubuntu/.ssh/instance_1.pem -r $1 ubuntu@146.118.64.185:$2

