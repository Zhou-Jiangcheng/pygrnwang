import pickle

from mpi4py import MPI
from multiprocessing import Pool

from pygrnwang.create_spgrn import *
from pygrnwang.utils import group


def pre_process(
    processes_num,
    event_depth_list,
    path_green,
    path_bin,
    spec_time_window,
    time_window,
    sampling_interval,
    before_tp,
    dist_range,
    delta_dist_range,
    max_frequency,
    path_nd=None,
    earth_model_layer_num=None,
):
    for event_depth in event_depth_list:
        create_dir(event_depth, path_green)
        create_inp(
            event_depth,
            path_green,
            spec_time_window,
            time_window,
            sampling_interval,
            before_tp,
            dist_range,
            delta_dist_range,
            max_frequency,
            path_nd,
            earth_model_layer_num,
        )
        path_bin_call = os.path.join(path_green, "spgrn2020.bin")
        if not os.path.exists(path_bin_call):
            shutil.copy(path_bin, path_bin_call)
    group_list = group(event_depth_list, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)
    return group_list


def create_grnlib_parallel_single_node(path_green):
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in group_list:
        print("computing " + str(item) + " km")
        for i in range(len(item)):
            item[i] = [item[i]] + [path_green]
        pool = Pool()
        r = pool.starmap_async(call_spgrn, item)
        r.get()
        pool.close()
        pool.join()
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def create_grnlib_parallel_multi_nodes(path_green):
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for ind_group in range(len(group_list)):
        comm = MPI.COMM_WORLD
        processes_num = comm.Get_size()
        rank = comm.Get_rank()
        if processes_num != len(group_list[0]):
            raise ValueError(
                "processes_num is %d, item num in group is %d. \n"
                "Pleasse check the process num!" % (
                    processes_num, len(group_list[0]))
            )
        print("ind_group:%d rank:%d" % (ind_group, rank))
        call_spgrn(event_depth=group_list[ind_group]
                   [rank], path_green=path_green)


if __name__ == "__main__":
    path_green_ = "/home/zjc/spgrn2020_lib/ak135fc_0_2"
    path_bin_ = "/home/zjc/python_works/pygrnwang/spgrn2020.bin"

    pre_process(
        processes_num=4,
        event_depth_list=[10 * h for h in range(8, 70)],
        path_green=path_green_,
        path_bin=path_bin_,
        spec_time_window=2048,
        time_window=2048,
        sampling_interval=1,
        before_tp=20,
        dist_range=[30, 100],
        delta_dist_range=[10, 10],
        max_frequency=0.2,
        path_nd="/home/zjc/python_works/pygrnwang/ak135.nd",
        earth_model_layer_num=34,
    )
    create_grnlib_parallel_single_node(path_green=path_green_)
    # create_grnlib_parallel_multi_nodes(path_green=path_green_)
