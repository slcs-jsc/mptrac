#! /bin/bash

# Create plots...
dot -Tpng clusters.dot -o clusters.png
dot -Tpdf clusters.dot -o clusters.pdf
