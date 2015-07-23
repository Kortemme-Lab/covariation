import sys
import os

metric_files = {
    'covariation_similarity' : 'covariation_similarity_backrub.txt',
    'profile_similarity' : 'backrub_sequence_profile.txt',
    'sequence_recovery' : 'sequence_recovery_backrub.txt',
    'sequence_entropy' : 'sequence_entropy_backrub.txt'
}
metric_R_names = {
    'covariation_similarity' : 'covariation_similarity',
    'profile_similarity' : 'profile_similarity',
    'sequence_recovery' : 'sequence_recovery',
    'sequence_entropy' : 'sequence_entropy'
}

this_dir = os.path.dirname(os.path.abspath(__file__))
for v in metric_files.values():
    assert(os.path.exists(os.path.join(this_dir, v)))

domains = ['PF00013','PF00018','PF00041','PF00072','PF00076','PF00085','PF00111','PF00168','PF00169','PF00179','PF00226','PF00240','PF00249','PF00254','PF00313','PF00327','PF00355','PF00364','PF00381','PF00439','PF00486','PF00498','PF00542','PF00550','PF00581','PF00582','PF00595','PF00691','PF00708','PF01029','PF01035','PF01451','PF01627','PF01833','PF02823','PF04002','PF07679','PF07686','PF08666','PF12844']
assert(len(domains) == 40)

published_methods = ['Fixed', '0.3', '0.6', '0.9', '1.2', '1.8', '2.4']


def get_method_ids(methods, backrub_method_prefix = '', backrub_method_suffix = ''):
    published_method_names = {'Fixed' : 'Fixed'}
    for m in published_methods[1:]:
        published_method_names[m] = backrub_method_prefix + m + backrub_method_suffix
    return [published_method_names[m] for m in methods]


def get_data_frame(methods, metric, backrub_method_prefix = '', backrub_method_suffix = ''):
    metric = metric.lower()

    if metric not in metric_files:
        raise Exception('The metric "{0}" is invalid. Please use one of "{1}".'.format(metric, '", "'.join(metric_files.keys())))

    all_metric_values = None
    try:
        metric_file = None
        metric_file = os.path.join(this_dir, metric_files[metric])
        all_metric_values = open(metric_file).read().strip()
    except:
        raise Exception('An error occurred trying to read the metric file "{0}".'.format(metric_file))

    metric_value_lines = [l.strip() for l in all_metric_values.split('\n') if l.strip()]
    assert(metric_value_lines[0].split() == ['Domain'] + published_methods)

    method_data = {}
    for method in methods:
        method_data[method] = []
        index = published_methods.index(method)
        for l in metric_value_lines[1:]:
            tokens = l.split()
            assert(len(tokens) == len(published_methods) + 1)
            method_data[method].append(tokens[index + 1])
        assert(len(method_data[method]) == 40)
        method_data[method] = ','.join(map(str, method_data[method]))

    for m in methods:
        if m not in published_methods:
            raise Exception('The results for method "{0}" are not available here. Please use one of "{1}".'.format(m, '", "'.join(published_methods)))

    published_method_names = {'Fixed' : 'Fixed'}
    for m in published_methods[1:]:
        published_method_names[m] = backrub_method_prefix + m + backrub_method_suffix

    num_methods = len(methods)

    all_domains_string = ",".join(["'{0}'".format("','".join(domains)) for m in methods])


    data_frame_string = '''
plos_benchmark_data <- data.frame(
  benchmark = c(''' + (','.join(["'Rosetta 3.2'" for x in range(num_methods * 40)])) + '''),
  method = c(''' + "'" + ("','".join(["','".join([published_method_names[m] for x in range(40)]) for m in methods])) + """'),
  domain = c(""" + all_domains_string + """),
  """ + metric_R_names[metric] + """ = c(\n""" + ',\n'.join([method_data[m] for m in methods]) + '''))'''

    return data_frame_string

