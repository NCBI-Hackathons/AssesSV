#!/bin/bash
samtools view -F 256 $1 | awk '{print length($10)}'
