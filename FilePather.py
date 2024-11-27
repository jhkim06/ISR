import json
from helper import get_bottom_dict


class FilePather:
    def __init__(self, base_dir, path_json):
        # TODO make a script to list samples and save as json file
        # Note that json file should in base_dir
        # TODO check if it is true!
        self.base_dir = base_dir
        self.file_path = dict()

        with open(self.base_dir + path_json, 'r') as config_file:
            config_json = json.load(config_file)
            self.file_path = config_json
        # self._update_file_path()

    def _update_file_path(self, temp_dict):
        # temp_dict is not nested dictionary
        for key in temp_dict:
            file_path = temp_dict[key]
            file_path = self.base_dir + file_path
            temp_dict[key] = file_path

    def get_path_dict(self,  *keys_tuple):

        out_dict = {}
        for keys in keys_tuple:
            file_path = self.file_path

            temp_dict = {}
            key_prefix = ''
            last_key = ''
            for key in keys:
                file_path = file_path[key]
                if key_prefix == '':
                    key_prefix = key
                else:
                    key_prefix += '_' + key
                last_key = key

            get_bottom_dict(file_path, temp_dict, key_prefix=key_prefix, last_key=last_key)
            self._update_file_path(temp_dict)
            out_dict.update(temp_dict)

        return out_dict
