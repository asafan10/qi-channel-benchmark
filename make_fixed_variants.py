"""Regenerate the corrected variants of our earlier implementation used for Table 5 (asp_fix.py):
qrc_fix1.py = quantum_ran_complete.py with the exact single-flip energy difference (Fault 5, Eq. (A3)) in DA and PA;
qrc_fix2.py = qrc_fix1.py with the standard replica-exchange criterion in DA (Fault 6, Eq. (18))."""
src = open("quantum_ran_complete.py").read()
old_da = "return float((1 - 2 * x[i]) * (col @ x - Q[i, i] * x[i]))"
old_pa = "(col @ pop[r] - Qs[i, i] * pop[r, i]))"
assert src.count(old_da) == 1 and src.count(old_pa) == 1
f1 = src.replace(old_da, "return float((1 - 2 * x[i]) * (col @ x - Q[i, i] * (2 * x[i] - 1)))") \
        .replace(old_pa, "(col @ pop[r] - Qs[i, i] * (2 * pop[r, i] - 1)))")
open("qrc_fix1.py", "w").write(f1)
old_sw = "dE   = energies[r + 1] - energies[r]"
assert f1.count(old_sw) == 1
open("qrc_fix2.py", "w").write(f1.replace(old_sw, "dE   = energies[r] - energies[r + 1]"))
print("wrote qrc_fix1.py, qrc_fix2.py")
