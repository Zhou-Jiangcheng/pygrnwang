import os
import pathlib

import obspy
import numpy as np
from obspy.taup import TauPyModel, taup_create
from .geo import d2km

def read_or_create_model(model_name, npz_dir=None):
    if model_name[-3:] == '.nd':
        path_nd = model_name
        model_name = os.path.splitext(os.path.basename(model_name))[0]
        try:
            model = TauPyModel(model=model_name)
        except:
            taup_create.build_taup_model(path_nd, output_folder=npz_dir, verbose=False)
            if npz_dir is not None:
                model_name = os.path.join(npz_dir, model_name + '.npz')
            model = TauPyModel(model=model_name)
    else:
        if npz_dir is not None:
            model_name = os.path.join(npz_dir, model_name + '.npz')
        try:
            model = TauPyModel(model=model_name)
        except:
            raise ValueError('You should provide nd file, or create model_name.npz file first!')
    return model

def remove_npz_file(model_name):
    if model_name[-3:] == '.nd':
        model_name = os.path.splitext(os.path.basename(model_name))[0]
    taup_data_dir = pathlib.Path(obspy.taup.__file__).parent / "data"
    npz_path = taup_data_dir / f"{model_name}.npz"
    os.remove(npz_path)

def cal_first_p(event_depth_km, dist_km, receiver_depth_km=0, model_name="ak135"):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    phases_list_p = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
    model = read_or_create_model(model_name)
    arrivals = model.get_travel_times(source_depth_in_km=event_depth_km,
                                      distance_in_degree=dist_km/d2km,
                                      phase_list=phases_list_p,
                                      receiver_depth_in_km=receiver_depth_km)
    return arrivals[0].time


def cal_first_p_s(event_depth_km, dist_km, receiver_depth_km=0, model_name="ak135"):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    phases_list_p = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
    model = read_or_create_model(model_name)
    arrivals_p = model.get_travel_times(source_depth_in_km=event_depth_km,
                                        distance_in_degree=dist_km/d2km,
                                        phase_list=phases_list_p,
                                        receiver_depth_in_km=receiver_depth_km)
    first_p = arrivals_p[0].time

    phases_list_s = ["s", "S", "sS", "pS", "Sg", "Sn", "Sdiff", "SKS"]
    arrivals_s = model.get_travel_times(source_depth_in_km=event_depth_km,
                                        distance_in_degree=dist_km/d2km,
                                        phase_list=phases_list_s,
                                        receiver_depth_in_km=receiver_depth_km)
    first_s = arrivals_s[0].time
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
