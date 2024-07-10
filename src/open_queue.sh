#!/bin/bash

srun --job-name=nb_compiler --pty --nodes=1  --cpus-per-task=8 -p gpu --gpus=1 -C a100 bash

