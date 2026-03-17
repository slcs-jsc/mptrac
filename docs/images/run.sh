#! /bin/bash

# Create plots...
dot -Tpng mptrac_modules.dot -o mptrac_modules.png
dot -Tpdf mptrac_modules.dot -o mptrac_modules.pdf
