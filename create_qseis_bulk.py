import pickle
from multiprocessing import Pool
from mpi4py import MPI

from pygrnwang.create_qseis import *
from pygrnwang.utils import group


def pre_process(
    processes_num,
    event_depth_list,
    path_green,
    time_window,
    sampling_interval,
    dist_range,
    delta_dist,
    isurf,
    rm_down=False,
    earth_model_layer_num=139,
    N_each_group=100,
    time_reduce_slowness=8,
):
    N_dist_group = None
    sub_sub_dirs_list = []
    for event_depth in event_depth_list:
        print("creating green func dir, info, inp for event_depth=%d" % event_depth)
        sub_dir, sub_sub_dirs, N_dist, N_dist_group = create_dir(
            event_depth, path_green, dist_range, delta_dist, N_each_group
        )
        create_greeninfo(
            event_depth,
            time_window,
            sampling_interval,
            dist_range,
            delta_dist,
            sub_dir,
            N_dist,
            N_dist_group,
            N_each_group,
        )
        create_inp(
            event_depth,
            path_green,
            time_window,
            sampling_interval,
            dist_range,
            delta_dist,
            sub_sub_dirs,
            N_dist,
            N_dist_group,
            isurf,
            rm_down,
            earth_model_layer_num,
            N_each_group,
            time_reduce_slowness,
        )
        sub_sub_dirs_list.append(sub_sub_dirs)
    jobs = []
    for i in range(len(event_depth_list)):
        for j in range(N_dist_group):
            jobs.append([event_depth_list[i], j])
    group_list = group(jobs, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)
    return group_list


def create_grnlib_parallel_single_node(path_green):
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in group_list:
        print("computing " + str(item))
        for i in range(len(item)):
            item[i] = item[i] + [path_green]
        pool = Pool()
        r = pool.starmap_async(call_qseis, item)
        r.get()
        pool.close()
        pool.join()

    e = datetime.datetime.now()
    print("total run time:" + str(e - s))


def create_grnlib_parallel_multi_nodes(path_green):
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    comm = MPI.COMM_WORLD
    processes_num = comm.Get_size()
    if processes_num != len(group_list[0]):
        raise ValueError(
            "processes_num is %d, item num in group is %d. \n"
            "Pleasse check the process num!" % (processes_num, len(group_list[0]))
        )
    rank = comm.Get_rank()
    ind_group = rank // processes_num
    ind_para = rank - ind_group * processes_num
    print(rank, ind_group, ind_para)
    call_qseis(
        event_depth=group_list[ind_group][ind_para][0],
        n_group=group_list[ind_group][ind_para][1],
        path_green=path_green,
    )


if __name__ == "__main__":
    pass
