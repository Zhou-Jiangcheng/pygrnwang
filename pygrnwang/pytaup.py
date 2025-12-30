import os
import sys

import numpy as np
import jpype.imports  # use jpype to call java class

if sys.platform == "win32":
    jar_path = os.path.join(sys.exec_prefix, 'Scripts', 'TauP.jar')
else:
    jar_path = os.path.join(sys.exec_prefix, 'bin', 'TauP.jar')

if not os.path.exists(jar_path):
    print(f"TauP.jar not found in {sys.exec_prefix}/bin or {sys.exec_prefix}/Scripts")
    print("Please install the TauP toolkit and ensure TauP.jar is in the correct directory.")
    sys.exit(1)

if not jpype.isJVMStarted():
    jpype.startJVM("--enable-native-access=ALL-UNNAMED", classpath=[jar_path])
from edu.sc.seis.TauP import TauP_Time  # type: ignore


def taup_time_java(
    event_depth_km, dist_km, phases_list, receiver_depth_km=0, model_name="ak135"
):
    ttobj = TauP_Time(model_name)
    ttobj.setSourceDepth(event_depth_km)
    ttobj.setReceiverDepth(receiver_depth_km)
    ttobj.setPhaseNames(phases_list)
    ttobj.calculate(dist_km / 111.19492664455874)

    N_arr = ttobj.getNumArrivals()
    results = {"phase": [], "puristphase": [], "time": [], "rayparameter": []}
    for i in range(N_arr):
        arr = ttobj.getArrival(i)
        results["phase"].append(str(arr.getName()))
        results["puristphase"].append(str(arr.getPuristName()))
        results["time"].append(float(arr.getTime()))
        results["rayparameter"].append(float(arr.getRayParam()))
    return results


def cal_first_p(event_depth_km, dist_km, receiver_depth_km=0, model_name="ak135"):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    phases_list_p = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
    ttobj_p = TauP_Time(model_name)
    ttobj_p.setSourceDepth(event_depth_km)
    ttobj_p.setReceiverDepth(receiver_depth_km)
    ttobj_p.setPhaseNames(phases_list_p)
    ttobj_p.calculate(dist_km / 111.19492664455874)
    arr_p = ttobj_p.getArrival(0)
    first_p = arr_p.getTime()
    return first_p


def cal_first_p_s(event_depth_km, dist_km, receiver_depth_km=0, model_name="ak135"):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    phases_list_p = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
    ttobj_p = TauP_Time(model_name)
    ttobj_p.setSourceDepth(event_depth_km)
    ttobj_p.setReceiverDepth(receiver_depth_km)
    ttobj_p.setPhaseNames(phases_list_p)
    ttobj_p.calculate(dist_km / 111.19492664455874)
    arr_p = ttobj_p.getArrival(0)
    first_p = arr_p.getTime()

    phases_list_s = ["s", "S", "sS", "pS", "Sg", "Sn", "Sdiff", "SKS"]
    ttobj_s = TauP_Time(model_name)
    ttobj_s.setSourceDepth(event_depth_km)
    ttobj_s.setReceiverDepth(receiver_depth_km)
    ttobj_s.setPhaseNames(phases_list_s)
    ttobj_s.calculate(dist_km / 111.19492664455874)
    arr_s = ttobj_s.getArrival(0)
    first_s = arr_s.getTime()

    return first_p, first_s


def create_tpts_table(
    path_green,
    event_depth_km,
    receiver_depth_km,
    dist_range,
    delta_dist,
    model_name="ak135",
    check_finished=False,
):
    path_tp_table = os.path.join(
        path_green,
        "%.2f" % event_depth_km,
        "%.2f" % receiver_depth_km,
        "tp_table.bin",
    )
    path_ts_table = os.path.join(
        path_green,
        "%.2f" % event_depth_km,
        "%.2f" % receiver_depth_km,
        "ts_table.bin",
    )
    if (
        check_finished
        and os.path.exists(path_tp_table)
        and os.path.exists(path_ts_table)
    ):
        return
    dist_kms = np.linspace(
        dist_range[0],
        dist_range[1],
        round(np.ceil((dist_range[1] - dist_range[0]) / delta_dist)) + 1,
    )
    tp_table = np.zeros(len(dist_kms), dtype=np.float32)
    ts_table = np.zeros(len(dist_kms), dtype=np.float32)
    for i in range(len(dist_kms)):
        first_p, first_s = cal_first_p_s(
            event_depth_km, dist_kms[i], receiver_depth_km, model_name
        )
        tp_table[i] = first_p
        ts_table[i] = first_s
    tp_table.tofile(path_tp_table)
    ts_table.tofile(path_ts_table)


def read_tpts_table(path_green, event_depth_km, receiver_depth_km, ind):
    fr_tp = open(
        os.path.join(
            path_green,
            "%.2f" % event_depth_km,
            "%.2f" % receiver_depth_km,
            "tp_table.bin",
        ),
        "rb",
    )
    tp = np.fromfile(file=fr_tp, dtype=np.float32, count=1, offset=ind * 4)[0]
    fr_tp.close()

    fr_ts = open(
        os.path.join(
            path_green,
            "%.2f" % event_depth_km,
            "%.2f" % receiver_depth_km,
            "ts_table.bin",
        ),
        "rb",
    )
    ts = np.fromfile(file=fr_ts, dtype=np.float32, count=1, offset=ind * 4)[0]
    fr_ts.close()
    return float(tp), float(ts)


if __name__ == "__main__":
    pass
