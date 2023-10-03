#!/bin/bash

 for i in $(cat locations.txt); do wget --user pippelshare --password \"'p@f0p!sh'\" "$i"; done 
