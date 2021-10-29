#!/usr/bin/bash

tabix -R coverage.bed amostra-lbb.HaplotypeCaller.PASS.vcf.gz > amostra-lbb.HaplotypeCaller.PASS.InBED.vcf.gz

tabix -p amostra-lbb.HaplotypeCaller.PASS.InBED.vcf.gz

