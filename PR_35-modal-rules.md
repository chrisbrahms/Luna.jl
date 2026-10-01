# gpu/35-modal-rules: fixed versus adaptive transverse integral

Base: `gpu/34-cpu-vs-gpu`. Benchmark scripts, results and report only; no change to `src/`.

`benchmark/modal_rules/`: five multimode cases (8 HE₁ₘ strong and weak, a self-compressing
HE₂ₘ vortex, mixed azimuthal orders, a rectangular guide), serial CPU runs of the adaptive
and fixed transverse rules against converged references, an overnight queue with
deadlines and PID watchdogs, and `REPORT.md`.

Result: do not switch the default to the fixed rule everywhere. The fixed rule is as
accurate and faster for radially symmetric, θ-structured and Cartesian sets, but needs an
estimated 2.5–3× the adaptive rule's time for the same spectral accuracy in the vortex
case. The adaptive default's `mfcn = 512` cap binds silently in every 2-D case (0.46 dB off
in the rectangular case), which is the main defect found; recommendations in the report.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01JLyXeRJFXy3CczpvjZHJWW
