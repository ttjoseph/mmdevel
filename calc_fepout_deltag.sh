#!/bin/bash

vmd -dispdev text <<FOO
package require parsefep
parsefep -gc 0 -forward $1
quit
FOO
