# env.sh -- source this in EVERY new shell before building or running:
#   source env.sh
# Edit the two paths marked EDIT to match your machine.

# --- 1) ROOT ----------------------------------------------------------------
# Conda route (recommended): activate the env that has root + gmsh.
#   conda activate garfield
# Binary/source route: source ROOT's setup script instead, e.g.
#   source /opt/root/bin/thisroot.sh          # EDIT if you use this route

# If you used conda, uncomment:
# conda activate garfield

# --- 2) Garfield++ ----------------------------------------------------------
export GARFIELD_INSTALL="$HOME/garfieldpp/install"          # EDIT
if [ -f "$GARFIELD_INSTALL/share/Garfield/setupGarfield.sh" ]; then
  source "$GARFIELD_INSTALL/share/Garfield/setupGarfield.sh"
else
  echo "WARNING: setupGarfield.sh not found under $GARFIELD_INSTALL"
fi

# Elmer (from the PPA) and Gmsh (pip/conda/apt) are already on PATH once installed.

# --- quick check ------------------------------------------------------------
echo "root-config : $(command -v root-config || echo MISSING)"
echo "gmsh        : $(command -v gmsh        || echo MISSING)"
echo "ElmerGrid   : $(command -v ElmerGrid   || echo MISSING)"
echo "ElmerSolver : $(command -v ElmerSolver || echo MISSING)"
echo "Garfield    : ${GARFIELD_INSTALL:-unset}"
