"""Reference implementation of the grounded labeling of a claim graph and a
randomized check that incremental recomputation equals full recomputation.

Nodes are claims and inference steps. Dependency relation: premise -> step,
step -> conclusion, objector -> target.
"""
import random
import sys

IN, OUT, UND = "in", "out", "und"


def make_graph(rng, n_claims, n_steps, n_obj):
    claims = list(range(n_claims))
    base = {c for c in claims if rng.random() < 0.4}
    steps = {}
    for s in range(n_steps):
        k = rng.randint(1, min(3, n_claims))
        steps[("s", s)] = (frozenset(rng.sample(claims, k)), rng.choice(claims))
    nodes = claims + list(steps)
    objs = set()
    for _ in range(n_obj):
        x, y = rng.choice(nodes), rng.choice(nodes)
        objs.add((x, y))
    return {"claims": set(claims), "base": base, "steps": steps, "objs": objs}


def signature(g):
    """Defining data of each node: the data its label depends on directly."""
    sig = {}
    sup = {c: set() for c in g["claims"]}
    for sid, (prem, conc) in g["steps"].items():
        sup[conc].add(sid)
    objectors = {}
    for x, y in g["objs"]:
        objectors.setdefault(y, set()).add(x)
    for c in g["claims"]:
        sig[c] = ("c", c in g["base"], frozenset(sup[c]),
                  frozenset(objectors.get(c, ())))
    for sid, (prem, conc) in g["steps"].items():
        sig[sid] = ("s", prem, conc, frozenset(objectors.get(sid, ())))
    return sig


def phi_node(g, sig, L, v):
    kind = sig[v][0]
    if kind == "s":
        prem, objs = sig[v][1], sig[v][3]
        if any(L[p] == OUT for p in prem) or any(L[o] == IN for o in objs):
            return OUT
        if all(L[p] == IN for p in prem) and all(L[o] == OUT for o in objs):
            return IN
        return UND
    base, sup, objs = sig[v][1], sig[v][2], sig[v][3]
    if any(L[o] == IN for o in objs):
        return OUT
    if not base and sup and all(L[s] == OUT for s in sup):
        return OUT
    if (base or any(L[s] == IN for s in sup)) and all(L[o] == OUT for o in objs):
        return IN
    return UND


def iterate(g, sig, nodes, L):
    """Iterate phi restricted to `nodes`, others held fixed in L, to a fixpoint."""
    L = dict(L)
    for v in nodes:
        L[v] = UND
    while True:
        new = {v: phi_node(g, sig, L, v) for v in nodes}
        if all(new[v] == L[v] for v in nodes):
            return L
        L.update(new)


def full(g):
    sig = signature(g)
    nodes = list(sig)
    return iterate(g, sig, nodes, {v: UND for v in nodes}), sig


def successors(g, sig):
    succ = {v: set() for v in sig}
    for v, d in sig.items():
        if d[0] == "s":
            for p in d[1]:
                succ[p].add(v)
            succ[v].add(d[2])
        for o in d[3]:
            succ[o].add(v)
    return succ


def incremental(g_new, old_labels, old_sig):
    sig = signature(g_new)
    affected = {v for v in sig if old_sig.get(v) != sig[v]}
    succ = successors(g_new, sig)
    U, stack = set(affected), list(affected)
    while stack:
        v = stack.pop()
        for w in succ[v]:
            if w not in U:
                U.add(w)
                stack.append(w)
    L = {v: old_labels[v] for v in sig if v not in U}
    for v in U:
        L[v] = UND
    return iterate(g_new, sig, list(U), L), U, sig


def mutate(rng, g):
    g2 = {"claims": set(g["claims"]), "base": set(g["base"]),
          "steps": dict(g["steps"]), "objs": set(g["objs"])}
    claims = sorted(g2["claims"])
    op = rng.choice(["add_obj", "del_obj", "add_step", "del_step",
                     "flip_base", "add_claim"])
    nodes = claims + list(g2["steps"])
    if op == "add_obj":
        g2["objs"].add((rng.choice(nodes), rng.choice(nodes)))
    elif op == "del_obj" and g2["objs"]:
        g2["objs"].discard(rng.choice(sorted(g2["objs"], key=str)))
    elif op == "add_step":
        k = rng.randint(1, min(3, len(claims)))
        sid = ("s", max([s[1] for s in g2["steps"]] + [-1]) + 1)
        g2["steps"][sid] = (frozenset(rng.sample(claims, k)), rng.choice(claims))
    elif op == "del_step" and g2["steps"]:
        sid = rng.choice(sorted(g2["steps"]))
        del g2["steps"][sid]
        g2["objs"] = {(x, y) for x, y in g2["objs"] if x != sid and y != sid}
    elif op == "flip_base":
        c = rng.choice(claims)
        g2["base"] ^= {c}
    elif op == "add_claim":
        g2["claims"].add(max(claims) + 1)
    return g2, op


def main(seed=0, trials=20000, updates=6):
    rng = random.Random(seed)
    checked = 0
    ops = {}
    cyc_like = 0
    for _ in range(trials):
        g = make_graph(rng, rng.randint(2, 8), rng.randint(0, 8), rng.randint(0, 8))
        labels, sig = full(g)
        for _ in range(updates):
            g2, op = mutate(rng, g)
            inc, U, sig2 = incremental(g2, labels, sig)
            ref, _ = full(g2)
            if inc != ref:
                print("MISMATCH", op, g2)
                sys.exit(1)
            checked += 1
            ops[op] = ops.get(op, 0) + 1
            g, labels, sig = g2, ref, sig2
    print("trials", trials, "updates checked", checked, "mismatches 0")
    print("operations", dict(sorted(ops.items())))


if __name__ == "__main__":
    main()
