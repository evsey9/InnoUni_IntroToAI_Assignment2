import mido
import music21
import random
import collections

# Constants
# Chords
chord_base = [0, 3, 7]
chord_major = chord_base.copy()
chord_minor = chord_base.copy()
chord_dim = chord_base.copy()
chord_sus2 = chord_base.copy()
chord_sus4 = chord_base.copy()
chord_major[1] = 4  # Assign chords like this to highlight changes from the minor choord, which is base
chord_dim[2] = 6
chord_sus2[1] = 2
chord_sus4[1] = 5
chord_list = [chord_major, chord_minor, chord_dim]
# Other
highest_note = 120


def get_midi_tempo(midi: mido.MidiFile) -> int:
    """
    Get the midi files' tempo
    :param midi: Midi file
    :return: Microseconds per beat (int)
    """
    for track in midi.tracks:
        for msg in track:
            if msg.type == "set_tempo":
                return msg.tempo
    return 500000  # Default tempo if no set_tempo messages


def get_velocity_mean(midi: mido.MidiFile) -> int:
    """
    Get the average velocity of a midi file
    :param midi: Midi file
    :return: Velocity (int)
    """
    count = 0
    velocity = 0
    for cur_note in midi.tracks[1]:  # only track 1 has notes
        if isinstance(cur_note, mido.messages.Message) and cur_note.type == "note_on":
            velocity += cur_note.velocity
            count += 1
    return velocity // count


def get_octave_mean(midi: mido.MidiFile) -> int:
    """
    Get the average octave of a midi file
    :param midi: Midi file
    :return: Octave (int)
    """
    count = 0
    octave = 0
    for cur_note in midi.tracks[1]:  # only track 1 has notes
        if isinstance(cur_note, mido.messages.Message) and cur_note.type == "note_on":
            octave += int(cur_note.note / 12)
            count += 1
    return octave // count


def get_note_length(track, chord_length):
    """
    Get total note length of a song
    :param track: The track of the song
    :param chord_length: Length of a chord
    :return:
    """
    time = 0
    for message in track:
        if type(message) is mido.messages.Message:
            time += message.time
    return (time + chord_length - 1) // chord_length


def get_notes(track, chord_length):
    """
    Get track notes
    :param track: The track with the messages
    :param chord_length: The length of a chord
    :return: All the notes of a track to use in chord computation
    """
    notes = [None for i in range(get_note_length(track, chord_length))]
    beats = 0
    final_note = 0
    for message in track:
        if type(message) is mido.messages.Message:
            cur_beat_index = beats // chord_length
            if message.type == "note_on":
                if beats % chord_length == 0 and notes[cur_beat_index] is None:
                    notes[cur_beat_index] = message.note
            elif message.type == "note_off":
                final_note = message.note
            beats += message.time
    if beats % chord_length == 0:
        notes[len(notes) - 1] = final_note
    return notes


class Chord:
    def __init__(self, base_note, c_type):
        """
        Chord class
        :param base_note: The note off which to build the chord
        :param c_type: The chord type (array of note displacements)
        """
        self.base_note = base_note
        self.chord_type = c_type
        self.notes = self.construct_list_from_type_and_note(base_note, c_type)

    @staticmethod
    def construct_list_from_type_and_note(base_note, c_type):
        """
        Constructs a list of notes in a chord from the base note and a chord type
        :param base_note: The note from which to construct chord
        :param c_type: The chord type (array of note displacements)
        :return: List of notes
        """
        return [(base_note + i_note) % 12 for i_note in c_type]

    def __contains__(self, item):
        """
        Returns if the chord contains a given note
        :param item: The note to check
        :return: True if the note is in the chord, false otherwise
        """
        for cur_note in self.notes:
            if cur_note == item % 12:
                return True
        else:
            return False

    def __eq__(self, other):
        """
        Check if chord is equal to another chord
        :param other: Other chord, or some other object
        :return: True if both chord type and base note are equal
        """
        if isinstance(other, Chord):
            return self.base_note == other.base_note and self.chord_type == other.chord_type
        else:
            return False  # A chord does not equal to a non-chord


class Accompaniment:
    def __init__(self, song_scale, song_size):
        """
        Accompaniment class, consisting of consonant chords, song scale and number of beats in the song
        :param song_scale: Song scale for determining consonant chords
        :param song_size: Number of beats in the song
        """
        self.song_scale = song_scale
        self.song_size = song_size
        self.consonant_chords = []
        self.generate_consonant_chords()

    def generate_consonant_chords(self):
        """
        Fills the consonant chords list with consonant chords
        """
        # Derive consonant chords from the circle of fifths:
        # Major chords:
        for i in range(5, 10, 2):
            self.consonant_chords.append(Chord((self.song_scale + (i % 9)) % 12, chord_major))
        # Minor chords:
        for i in range(7, 12, 2):
            self.consonant_chords.append(Chord((self.song_scale + (i % 9 + 2)) % 12, chord_minor))
        # Diminished chord:
        self.consonant_chords.append(Chord((self.song_scale + 11) % 12, chord_dim))

    def consonant_chords_contain_chord(self, chord):
        """
        Checks if the chord is in the accompaniment's consonant chords
        :param chord: The chord to check
        :return: True if the consonant chords list contains the chord, false otherwise
        """
        return any([consonant_chord_i == chord for consonant_chord_i in self.consonant_chords])

    def consonant_chord_notes_contain_note(self, note):
        """
        Checks if the note is contained in any of the consonant chords
        :param note: The note to check
        :return: True if the consonant chords contain the note
        """
        return any([note in consonant_chord_i for consonant_chord_i in self.consonant_chords])


class Chromosome:
    def __init__(self, chromosome_size):
        """
        Chromosome for genetic algorithm, representing an accompaniment of chords
        :param chromosome_size: Number of genes in chromosome
        """
        self.chromosome_size = chromosome_size
        self.gene_pool = [None for i in range(self.chromosome_size)]
        self.rating = 0
        self.generate_genes()

    def generate_genes(self):
        """
        Generate a set of random chord genes
        """
        self.gene_pool = [generate_random_chord() for i in range(self.chromosome_size)]

    def __eq__(self, other):
        """
        Compare by rating
        :param other: Other chromosome
        :return: True if both chromosomes have equal rating, false otherwise or if other object is not a chromosome
        """
        if isinstance(other, Chromosome):
            return self.rating == other.rating
        else:
            return False

    def __lt__(self, other):
        """
        Compare by rating
        :param other: Other chromosome
        :return: True if other chromosome has lesser rating, false otherwise
        """
        return self.rating < other.rating


def generate_random_note(lowest_note_local: int = 0, highest_note_local: int = highest_note):
    """
    Returns a random note from lowest_note_local to highest_note_local
    :param lowest_note_local: Lowest note possible
    :param highest_note_local: Highest note possible
    :return: Random note int
    """
    return random.randint(lowest_note_local, highest_note_local)


def generate_random_chord(lowest_note_local: int = 0, highest_note_local: int = highest_note):
    """
    Returns a random chord with base note from lowest_note_local to highest_note_local
    :param lowest_note_local: Lowest note possible
    :param highest_note_local: Highest note possible
    :return:
    """
    return Chord(generate_random_note(lowest_note_local, highest_note_local), random.choice(chord_list))


def generate_population(population_size, chromo_size):
    """
    Generate a population of size population_size with chromosome of chromo_size
    :param population_size:
    :param chromo_size:
    :return: Population
    """
    return [Chromosome(chromo_size) for i in range(population_size)]


def calculate_new_chromosome_ratings(chromosome, accompaniment: Accompaniment, notes):
    """
    Calculates the new rating for a chromosome
    :param chromosome: The chromosome itself
    :param accompaniment: Accompaniment with consonant chords
    :param notes: The song notes
    :return: New rating for the chromosome
    """
    new_rating = chromosome.chromosome_size
    for i, gene in enumerate(chromosome.gene_pool):
        if accompaniment.consonant_chords_contain_chord(gene):
            new_rating -= 0.5
            if notes[i] is None:
                new_rating -= 0.5
                continue
            if not accompaniment.consonant_chord_notes_contain_note(notes[i]) or notes[i] in gene:
                new_rating -= 0.5
    return new_rating


def adjust_ratings(population, accompaniment: Accompaniment, notes):
    """
    Adjusts the ratings of each chromosome in the population - the fitness function of the genetic algorithm
    :param notes: Song notes
    :param accompaniment: Accompaniment with consonant chords
    :param population: The population with the chromosomes
    """
    chromosome: Chromosome
    for chromosome in population:
        chromosome.rating = calculate_new_chromosome_ratings(chromosome, accompaniment, notes)


def select_top(population, survivor_list):
    """
    Selects the top chromosomes from the population
    :param population: The population
    :param survivor_list: The list of the survivors
    :return: New list of survivors
    """
    survivors_size = len(survivor_list)
    survivors = [population[i] for i in range(survivors_size)]
    return survivors


def get_random_chromosome_index(chromosomes, first_parent_index=None):
    """
    Gets a random chromosome index
    :param chromosomes: The list of chromosomes, used to know the size
    :param first_parent_index: The first parent, so that we don't take an already used index
    :return: Random index
    """
    chromosome_number = len(chromosomes)
    while True:
        rand_index = random.randrange(0, chromosome_number)
        if first_parent_index is None or first_parent_index != rand_index:
            return rand_index


def cross_chromosomes(first_chromosome: Chromosome, second_chromosome: Chromosome):
    """
    Slice the parent chromosomes at a random point and cross them over
    :param first_chromosome: First chromosome
    :param second_chromosome: Second chromosome
    :return: New chromosome
    """
    chromo_size = first_chromosome.chromosome_size
    slicing_point = random.randrange(0, chromo_size)
    child_chromosome = Chromosome(chromo_size)
    # Slice using list slices
    child_chromosome.gene_pool[:slicing_point] = first_chromosome.gene_pool[:slicing_point]
    child_chromosome.gene_pool[slicing_point:chromo_size] = second_chromosome.gene_pool[slicing_point:chromo_size]
    return child_chromosome


def repopulate(population, survivors, child_count):
    """
    Repopulate by crossing the survivors
    :param population: Population to fill
    :param survivors: The parents
    :param child_count: Current child count
    :return: New population
    """
    pop_size = len(population)
    while child_count < pop_size:
        first_parent_index = get_random_chromosome_index(survivors)
        second_parent_index = get_random_chromosome_index(survivors, first_parent_index)
        first_parent = survivors[first_parent_index]
        second_parent = survivors[second_parent_index]
        population[child_count] = cross_chromosomes(first_parent, second_parent)
        population[child_count + 1] = cross_chromosomes(second_parent, first_parent)
        child_count += 2
    return population


def mutate_chromosomes(population, chromosome_mutation_number, gene_mutation_number):
    """
    Mutates some genes in some chromosomes
    :param population: The population
    :param chromosome_mutation_number: How many chromosomes to mutate
    :param gene_mutation_number: How many genes in each chromosome to mutate
    :return: New population
    """
    population_size = len(population)
    for i in range(chromosome_mutation_number):
        chromosome_index = get_random_chromosome_index(population)
        chromosome = population[chromosome_index]
        for j in range(gene_mutation_number):
            new_chord = generate_random_chord()
            gene_index = random.randrange(0, chromosome.chromosome_size)
            chromosome.gene_pool[gene_index] = new_chord
    return population


def main():
    input_file = "input1.mid"
    output_file = "output1.mid"
    midi_file = mido.MidiFile(input_file, clip=True)
    midi_file.type = 1
    midi_score = music21.converter.parse(input_file)
    key = midi_score.analyze("key")  # use music21 to get the key
    print("Key is:", key)
    scale = (key.tonic.midi + (key.mode == "minor") * 3) % 12  # For consonant chord calculation
    chord_length = midi_file.ticks_per_beat * 2  # Two beats
    tempo = get_midi_tempo(midi_file)
    notes = get_notes(midi_file.tracks[1], chord_length)
    notes_length = len(notes)
    accompaniment_velocity = int(get_velocity_mean(midi_file) * 0.9)  # The volume of our accompaniment
    accompaniment_displacement = 12 * (get_octave_mean(midi_file) - 1)  # The range of our accompaniment
    accompaniment_genome = Accompaniment(scale, notes_length)

    track = [mido.MetaMessage("set_tempo", tempo=tempo, time=0)]  # The track that will contain our chords

    # Genetic algorithm:
    pop_size = 128  # Population size
    survivors = [None for i in range(pop_size // 4)]  # Best chromosomes
    population = generate_population(pop_size, accompaniment_genome.song_size)
    iteration_i = 0
    iteration_max = 5000
    # Select best of population, repopulate and mutate
    while True:
        adjust_ratings(population, accompaniment_genome, notes)
        population.sort()
        if population[0].rating == 0 or iteration_i >= iteration_max:
            break  # We either reached the perfect accompaniment or we have iterated for long enough
        survivors = select_top(population, survivors)
        population = repopulate(population, survivors, len(population) // 2)
        population = mutate_chromosomes(population, len(population) // 2, 1)

    # Writing
    chord: Chord
    for chord in population[0].gene_pool:
        for i in range(3):
            track.append(
                mido.messages.Message('note_on', channel=0, note=chord.notes[i] + accompaniment_displacement,
                                      velocity=accompaniment_velocity,
                                      time=0))
        for i in range(3):
            track.append(
                mido.messages.Message('note_off', channel=0, note=chord.notes[i] + accompaniment_displacement,
                                      velocity=accompaniment_velocity,
                                      time=(chord_length if i == 0 else 0)))
    midi_file.tracks.append(track)
    midi_file.save(output_file)


if __name__ == "__main__":
    main()
