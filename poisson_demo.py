"""Chạy ví dụ Poisson chính thức của scikit-fem và kiểm tra nghiệm.

Nguồn: https://scikit-fem.readthedocs.io/en/stable/gettingstarted.html
Miền: (0, 1) x (0, 1); -Laplace(u) = sin(pi*x) sin(pi*y); u = 0 trên biên.
Phần giải bám theo hướng dẫn; bổ sung lưu hình, dữ liệu và kiểm tra lưới.
Chạy: python poisson_demo.py
"""

from pathlib import Path
import csv
import json
import platform

import numpy as np
import scipy
import skfem as fem
from skfem.helpers import dot, grad
import matplotlib

# Cho phép lưu ảnh khi chạy trên cloud không có cửa sổ đồ họa.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


@fem.BilinearForm
def stiffness(u, v, _):
    """Vế trái dạng yếu: tích vô hướng hai gradient."""
    return dot(grad(u), grad(v))


@fem.LinearForm
def load(v, w):
    """Vế phải dạng yếu: f(x,y) nhân hàm thử v."""
    x, y = w.x
    return np.sin(np.pi * x) * np.sin(np.pi * y) * v


def exact_solution(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y) / (2.0 * np.pi**2)


@fem.Functional
def squared_error(w):
    x, y = w.x
    return (w["uh"] - exact_solution(x, y)) ** 2


def solve_poisson(refinement):
    # 1. Tạo lưới tam giác của hình vuông đơn vị.
    mesh = fem.MeshTri().refined(refinement)

    # 2. Dùng phần tử tam giác P1: hàm tuyến tính trên từng tam giác.
    basis = fem.Basis(mesh, fem.ElementTriP1())

    # 3. Lắp ráp hệ phương trình A * solution = b.
    A = stiffness.assemble(basis)
    b = load.assemble(basis)

    # 4. Đặt u = 0 trên toàn bộ biên, loại ẩn biên và giải hệ.
    boundary = basis.get_dofs()
    solution = fem.solve(*fem.condense(A, b, D=boundary))

    # 5. Kiểm tra với nghiệm chính xác trong tài liệu.
    error2 = float(squared_error.assemble(basis, uh=basis.interpolate(solution)))
    l2 = float(np.sqrt(error2))
    # Chuẩn L2 chính xác của nghiệm trên hình vuông bằng 1/(4*pi^2).
    relative_l2 = l2 / (1.0 / (4.0 * np.pi**2))
    free = np.setdiff1d(np.arange(basis.N), boundary.all())
    residual = (A @ solution - b)[free]
    rel_residual = float(np.linalg.norm(residual) / np.linalg.norm(b[free]))
    stats = {
        "refinement_level": refinement,
        "triangles": int(mesh.nelements),
        "nodes": int(mesh.nvertices),
        "degrees_of_freedom": int(basis.N),
        "boundary_nodes": int(len(boundary.all())),
        "integrated_squared_error": error2,
        "l2_error": l2,
        "relative_l2_error_percent": 100.0 * relative_l2,
        "max_boundary_value_abs": float(np.max(np.abs(solution[boundary.all()]))),
        "relative_linear_system_residual": rel_residual,
        "max_computed_value": float(solution.max()),
        "max_exact_value": float(1.0 / (2.0 * np.pi**2)),
    }
    return mesh, solution, stats


def main():
    folder = Path(__file__).resolve().parent
    # Lưới mức 3 đúng ví dụ Getting started; mức 4 và 5 để kiểm tra hội tụ.
    runs = [solve_poisson(level) for level in (3, 4, 5)]
    mesh, solution, main_stats = runs[0]
    all_stats = [item[2] for item in runs]

    # Các kiểm tra này đối chiếu điều kiện toán học và số liệu hướng dẫn.
    assert main_stats["triangles"] == 128
    assert main_stats["nodes"] == 81
    assert np.isclose(main_stats["integrated_squared_error"], 1.069e-6, rtol=0.001)
    for stats in all_stats:
        assert stats["max_boundary_value_abs"] < 1e-12
        assert stats["relative_linear_system_residual"] < 1e-10
    for coarse, fine in zip(all_stats, all_stats[1:]):
        assert fine["l2_error"] < coarse["l2_error"]
        fine["l2_reduction_factor"] = coarse["l2_error"] / fine["l2_error"]

    # Lưu hình lưới và nghiệm của trường hợp chính (128 tam giác).
    tri = mtri.Triangulation(mesh.p[0], mesh.p[1], mesh.t.T)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), constrained_layout=True)
    fig.suptitle("Ví dụ Poisson 2D bằng scikit-fem", fontsize=16, fontweight="bold")
    axes[0].triplot(tri, color="#53647A", linewidth=0.65)
    axes[0].set_title("Lưới: 128 tam giác, 81 nút", fontsize=12)
    filled = axes[1].tripcolor(tri, solution, shading="gouraud", cmap="viridis")
    axes[1].set_title("Nghiệm tính được — u = 0 trên biên", fontsize=12)
    colorbar = fig.colorbar(filled, ax=axes[1], shrink=0.82, pad=0.03)
    colorbar.set_label("u(x, y)")
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    fig.savefig(folder / "Ket_qua_Poisson.png", dpi=200, facecolor="white")
    plt.close(fig)

    # Dữ liệu tại nút, thuận tiện mở bằng bảng tính hoặc xử lý tiếp.
    exact_nodes = exact_solution(mesh.p[0], mesh.p[1])
    np.savetxt(
        folder / "nodal_solution.csv",
        np.column_stack([mesh.p[0], mesh.p[1], solution, exact_nodes, solution-exact_nodes]),
        delimiter=",", header="x,y,u_fem,u_exact,nodal_error", comments="",
    )
    keys = list(all_stats[-1])
    with (folder / "mesh_convergence.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(all_stats)
    results = {
        "source": "https://scikit-fem.readthedocs.io/en/stable/gettingstarted.html",
        "versions": {
            "python": platform.python_version(), "scikit-fem": fem.__version__,
            "numpy": np.__version__, "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "primary_case": main_stats,
        "mesh_convergence": all_stats,
        "checks_passed": True,
        "error_note": "The tutorial value 1.069e-6 is the integrated squared error; take its square root for the L2 norm.",
    }
    (folder / "results.json").write_text(json.dumps(results, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
