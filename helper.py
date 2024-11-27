# Note structure of dictionary represents directory structure
def get_bottom_dict(nested_dict, out_dict, key_history='', key_prefix='', last_key=''):
    # key_history: search history
    # key_prefix: starting point of searching

    if isinstance(nested_dict, dict):
        for new_key in nested_dict:
            if key_history == '':
                new_key_history = new_key
            else:
                new_key_history = key_history + '_' + new_key
            # recursive call
            get_bottom_dict(nested_dict[new_key], out_dict, new_key_history, key_prefix, new_key)
    # Reached bottom of nested dictionary
    else:
        path_from_key_history = ''
        if key_history[-len('_' + last_key):] == '_' + last_key:
            path_from_key_history = '/'.join(key_history[:-len('_' + last_key)].split('_')) + '/'

        # started from the top of nested dictionary
        if key_prefix == '':
            out_dict[key_history] = path_from_key_history + nested_dict
        else:
            path_from_key_prefix = key_prefix
            if key_history == '' and last_key != '':
                path_from_key_prefix = path_from_key_prefix[:-len('_'+last_key)]
            path_from_key_prefix = '/'.join(path_from_key_prefix.split('_')) + '/'

            if key_history == '':
                out_dict[key_prefix] = path_from_key_prefix + path_from_key_history + nested_dict
            else:
                out_dict[key_prefix + '_' + key_history] = path_from_key_prefix + path_from_key_history + nested_dict
