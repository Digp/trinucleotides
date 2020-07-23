#!/bin/bash

date
R CMD build minipack
R CMD INSTALL minipack*
date
