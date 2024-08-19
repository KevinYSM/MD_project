#!/bin/bash
sudo scp -i /home/ubuntu/.ssh/kevin_instance.pem -r $1 ubuntu@146.118.68.29:$2
