import random
from numpy import array
from numpy import sum


def bootstrap(data, sample_amount, size):
	samples = []
	for iteration in range(size):
		sample = []
		for i in range(sample_amount):
			sample.append(random.choice(data))
		sample = array(sample)
		temp = sum(sample, 0)
		samples.append(list(temp))
	return samples


def bootstrap_scores(samples, number_of_kmeres, number_of_filters):
	scores = []
	# calculates float for each value in [hits per filter]
	for i in range(len(samples)):
		score = []
		for j in range(number_of_filters):
			if samples[i][j] == 0:
				score.append(0.0)
			else:
				score.append(round(float(samples[i][j]) / float(number_of_kmeres), 2))
		scores.append(score)
	return scores

def cluster_kmers(kmer_list, kmer_dict):
    clusters = {}
    contig_median_list = []
    # Schleife über alle kmere in kmer_list
    for i in range(len(kmer_list)):
        kmer = kmer_list[i]
        # Überprüfen, ob das kmer in kmer_dict vorhanden ist
        kmer_info = kmer_dict.get(kmer)
        if kmer_info is None:
            continue
        # Holt die Contig-ID und kmer-Position aus kmer_dict
        kmer_id = kmer_dict[kmer][1]
        kmer_pos = kmer_dict[kmer][0]
        # Füge das kmer dem entsprechenden Contig in das Dictionary hinzu
        # Ein Contig ist ein cluster
        if kmer_id not in clusters:
            clusters[kmer_id] = []
        clusters[kmer_id].append((kmer, kmer_pos))
    # Schleife über alle Contigs im Dictionary clusters
    for contig in clusters:
        contig_list = clusters[contig]
        contig_len = len(contig_list)
        if contig_len < 2:
            #print("Zu wenig kmere im Contig!")
            continue
        # Sortieren der kmere in der Liste nach der Position
        sorted_contig_list = sorted(contig_list, key=lambda x: x[1])
        distances = []
        # Schleife über die sortierte Liste von kmere im Contig
        for i in range(1, len(sorted_contig_list)):
            kmer_pos = sorted_contig_list[i][1]
            prev_kmer_pos = sorted_contig_list[i-1][1]
            # Berechnung der Distanz zum vorherigen kmer
            distance = kmer_pos - prev_kmer_pos
            distances.append(distance)
        # Berechnung der Median-Entfernung für das aktuelle Contig
        #print(distances)
        median_distance = statistics.median(distances)
        #print(median_distance)
        contig_median_list.append((contig, median_distance, contig_len))
    # Summe der Contigs in contig_median_list
    num_contigs = len(contig_median_list)
    # Liste aller Contig-Größen
    contig_lengths = [x[2] for x in contig_median_list]
    # Liste aller Median-Entfernungen
    median_distances = [x[1] for x in contig_median_list]
    # Median der Median-Entfernungen berechnen
    if len(median_distances) > 0:
        median_of_medians = statistics.median(median_distances)
    else:
        median_of_medians = None
    # Ergebnisliste erstellen
    result = [num_contigs, median_of_medians, contig_lengths]
    return result