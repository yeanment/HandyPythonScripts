#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
"""
Implemented to remove calculated OpenFOAM case file after certain time. 
"""

import sys, os
import shutil

def removeAfterTime1(time, path_root):
    # os.chdir(path_root)
    for dir_path, dir_names, file_names in os.walk(path_root):
        for dir_name in dir_names:
            try:
                dir_time = float(dir_name)
                if (dir_time > time):
                    path_remove = os.path.join(dir_path, dir_name)
                    print("Remove dir {0}.".format(path_remove))
                    shutil.rmtree(path_remove)
            except ValueError:
                pass
            else:
                pass
    return

def removeAfterTime(time, path_root, maxdepth = None):
    try:
        names = os.listdir(path_root)
        names.sort()
    except OSError as e:
        print("OSError in removeAfterTime: ", e)
        return
    else:
        pass
    dirs, files = [], []
    for name in names:
        if os.path.isdir(os.path.join(path_root, name)):
            dirs.append(name)
            try:
                dir_time = float(name)
                if (dir_time > time):
                    path_remove = os.path.join(path_root, name)
                    print("Remove dir {0}.".format(path_remove))
                    shutil.rmtree(path_remove)
                    dirs.pop()
            except ValueError as e:
                # print("OSError in removeAfterTime: ", e)
                pass
            else:
                pass
        else:
            files.append(name)
    
    yield path_root, dirs, files
    
    if maxdepth is None or maxdepth > 1:
        for dir_name in dirs:
            dir_path = os.path.join(path_root, dir_name)
            for x in removeAfterTime(time, dir_path, None if maxdepth is None else maxdepth-1):
                yield x


def contained_dirs(dir):
    return filter(os.path.isdir,
                  [os.path.join(dir, f) for f in os.listdir(dir)])

if __name__ == "__main__":
    if '--time' in sys.argv:
        time_string_index = sys.argv.index('--time') + 1
        try:
            time_string = sys.argv[time_string_index]
            time_after = float(time_string)
            print(time_after)
        except IndexError:
            time_after = 0.01
        except ValueError:
            print('Error in input time.')
            exit()
        else:
            pass
    else:
        time_after = 0.01
        print('************')
        print('Remove calculated file after 0.01.')
        print('************')
    # print('ok')
    if '--case' in sys.argv:
        path_root_index = sys.argv.index('--case') + 1
        try:
            path_root = sys.argv[path_root_index]    
        except IndexError:
            path_root = '.'
        else:
            pass
    else:
        path_root = '.'
        print('************')
        print('Remove calculated file in .')
        print('************')
    
    if '--depth' in sys.argv:
        depth_string_index = sys.argv.index('--depth') + 1
        try:
            depth_string = sys.argv[depth_string_index]
            depth = int(depth_string)
        except IndexError:
            depth = 1
        except ValueError:
            print('Error in input depth.')
            exit()
        else:
            pass
    else:
        depth = 1
        print('************')
        print('Remove calculated file with depth 1.')
        print('************')
    
    path_root = os.path.abspath(path_root)
    if (os.access(path_root, os.F_OK) == False):
        print('Path not exists.')
    if os.path.isdir(path_root) == False:
        print('Not a path.')
    print('**********************')
    print("Remove case file after time {0} in path {1} with depth {2}.".format(time_after, path_root,depth))
    
    for path, dirs, files in removeAfterTime(time_after, path_root, depth):
        print("Executing in {0}".format(path))
