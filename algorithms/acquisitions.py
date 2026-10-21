"""Preset batch acquisition functions that a run can select by name.

The catalog is plain data so the web server can list and validate presets
without importing torch; BoTorch is imported only inside the builders. Every
preset is a Monte Carlo acquisition, so sequential greedy batch selection in
optimize_acqf_discrete (which sets X_pending) works for any batch size.
Models are fitted on direction-signed targets, so every preset maximizes.
"""
from __future__ import annotations

PARAMETERS = {
    "ucb_beta": {"label": "UCB β", "default": 2.0, "min": 0.0, "max": 100.0,
                 "description": "探索权重：越大越偏向不确定区域，0 等价于后验均值"},
}


def _sampler(args):
    import torch
    from botorch.sampling.normal import SobolQMCNormalSampler
    return SobolQMCNormalSampler(torch.Size([args.mc_samples]), seed=args.seed)


def _best_signed(data, args):
    return (data.train_Y.reshape(-1) * (-1 if args.minimize else 1)).max()


def _signs(data, args):
    return data.pool_X.new_tensor([-1 if d == "min" else 1 for d in args.directions])


def _qlognei(model, data, args, **_):
    from botorch.acquisition.logei import qLogNoisyExpectedImprovement
    return qLogNoisyExpectedImprovement(model=model, X_baseline=data.train_X,
                                        sampler=_sampler(args))


def _qlogei(model, data, args, **_):
    from botorch.acquisition.logei import qLogExpectedImprovement
    return qLogExpectedImprovement(model=model, best_f=_best_signed(data, args),
                                   sampler=_sampler(args))


def _qucb(model, data, args, *, ucb_beta, **_):
    from botorch.acquisition.monte_carlo import qUpperConfidenceBound
    return qUpperConfidenceBound(model=model, beta=ucb_beta, sampler=_sampler(args))


def _qpi(model, data, args, **_):
    from botorch.acquisition.monte_carlo import qProbabilityOfImprovement
    return qProbabilityOfImprovement(model=model, best_f=_best_signed(data, args),
                                     sampler=_sampler(args))


def _qsr(model, data, args, **_):
    from botorch.acquisition.monte_carlo import qSimpleRegret
    return qSimpleRegret(model=model, sampler=_sampler(args))


def _qstd(model, data, args, **_):
    from botorch.acquisition.monte_carlo import qPosteriorStandardDeviation
    return qPosteriorStandardDeviation(model=model, sampler=_sampler(args))


def _qlognehvi(model, data, args, ref_point, **_):
    from botorch.acquisition.multi_objective.logei import (
        qLogNoisyExpectedHypervolumeImprovement)
    return qLogNoisyExpectedHypervolumeImprovement(
        model=model, ref_point=data.pool_X.new_tensor(ref_point) * _signs(data, args),
        X_baseline=data.train_X, sampler=_sampler(args))


def _qlogehvi(model, data, args, ref_point, **_):
    from botorch.acquisition.multi_objective.logei import qLogExpectedHypervolumeImprovement
    from botorch.utils.multi_objective.box_decompositions.non_dominated import (
        FastNondominatedPartitioning)
    signs = _signs(data, args)
    ref = data.pool_X.new_tensor(ref_point) * signs
    partitioning = FastNondominatedPartitioning(ref_point=ref, Y=data.train_Y * signs)
    return qLogExpectedHypervolumeImprovement(model=model, ref_point=ref,
                                              partitioning=partitioning,
                                              sampler=_sampler(args))


def _qlognparego(model, data, args, **_):
    from botorch.acquisition.multi_objective.parego import qLogNParEGO
    return qLogNParEGO(model=model, X_baseline=data.train_X, sampler=_sampler(args))


ACQUISITIONS = {
    "single": {
        "qlognei": {"label": "qLogNEI", "params": [], "build": _qlognei,
                    "description": "噪声期望改进（默认）：以全部已测配方为 baseline，适合含噪声的实验"},
        "qlogei": {"label": "qLogEI", "params": [], "build": _qlogei,
                   "description": "期望改进：以已测最优值为基准，假定观测近似无噪声"},
        "qucb": {"label": "qUCB", "params": ["ucb_beta"], "build": _qucb,
                 "description": "置信上界：均值 + β 加权的不确定度，用 β 调节探索与利用"},
        "qpi": {"label": "qPI", "params": [], "build": _qpi,
                "description": "改进概率：只看超过已测最优的概率，偏保守、偏利用"},
        "qsr": {"label": "qSR", "params": [], "build": _qsr,
                "description": "后验均值（纯利用）：直接挑预测最好的配方"},
        "qstd": {"label": "qPosteriorStd", "params": [], "build": _qstd,
                 "description": "后验标准差（纯探索）：挑模型最不确定的配方，用于主动学习建模"},
    },
    "multi": {
        "qlognehvi": {"label": "qLogNEHVI", "params": [], "build": _qlognehvi,
                      "description": "噪声期望超体积改进（默认）：以全部已测配方为 baseline"},
        "qlogehvi": {"label": "qLogEHVI", "params": [], "build": _qlogehvi,
                     "description": "期望超体积改进：以已测 Pareto 前沿为基准，假定观测近似无噪声"},
        "qlognparego": {"label": "qLogNParEGO", "params": [], "build": _qlognparego,
                        "description": "随机权重 Chebyshev 标量化：目标较多时更快，不使用参考点"},
    },
}
DEFAULT_ACQUISITION = {"single": "qlognei", "multi": "qlognehvi"}


def catalog():
    """JSON-serializable presets for the web form."""
    return {
        "defaults": DEFAULT_ACQUISITION,
        "parameters": PARAMETERS,
        **{mode: [{"id": name, "label": spec["label"], "description": spec["description"],
                   "params": spec["params"]} for name, spec in presets.items()]
           for mode, presets in ACQUISITIONS.items()},
    }


def selected(args, mode):
    """Preset name and its parameter values; namespaces without the option get the default."""
    name = getattr(args, "acquisition", None) or DEFAULT_ACQUISITION[mode]
    params = {key: getattr(args, key, PARAMETERS[key]["default"])
              for key in ACQUISITIONS[mode][name]["params"]}
    return name, params


def build_acquisition(mode, model, data, args, ref_point=None):
    name, params = selected(args, mode)
    return ACQUISITIONS[mode][name]["build"](model, data, args, ref_point=ref_point, **params)
