import random
import pandas as pd
import numpy as np
from numpy.linalg import pinv
from numpy import matmul
from collections import Counter

from psychopy.visual import rect, circle, polygon
from psychopy import visual, logging, core

import sys
sys.path.append('.')
import walks
import graphs

keypress_map = {'1': 1,
                '2': 2,
                '3': 3,
                '4': 4,
                '5': 5,
                ' ': ' ',
                'j': 'j',
                'k': 'k',
                'l': 'l',
                ';': ';',
                'q': 'quit',
                'escape': 'quit'}

response_keys = [' ', 'j', 'k', 'l', ';']

FINGER_COMBINATIONS = [[True, False, False, False, False],
                       [False, True, False, False, False],
                       [False, False, True, False, False],
                       [False, False, False, True, False],
                       [False, False, False, False, True],
                       [True, True, False, False, False],
                       [True, False, True, False, False],
                       [True, False, False, True, False],
                       [True, False, False, False, True],
                       [False, True, True, False, False],
                       [False, True, False, True, False],
                       [False, True, False, False, True],
                       [False, False, True, True, False],
                       [False, False, True, False, True],
                       [False, False, False, True, True]]

target_offsets = [-0.1, 0.0, 0.0, 0.0, 0.0]

COLORS = dict(
    green=(27, 158, 119),  # green
    orange=(217, 95, 2),  # orange
    purple=(137, 132, 199),  # purple
    magenta=(231, 105, 138),  # magenta
    red=(208, 32, 32),  # red
    black=(0, 0, 0)
)

BG_COLOR = (100, 100, 100)

color_perms = [
    {'random': COLORS['orange'], 'modular': COLORS['green'], 'lattice': COLORS['magenta'], 'demo': COLORS['purple']},
    {'random': COLORS['orange'], 'modular': COLORS['magenta'], 'lattice': COLORS['green'], 'demo': COLORS['purple']}
]

cue_perms = [
    {'demo': 'square', 'random': 'triangle', 'modular': 'pentagon', 'lattice': 'circle'},
    {'demo': 'square', 'random': 'triangle', 'modular': 'circle', 'lattice': 'pentagon'}
]

class Target(rect.Rect):
    """

    """

    def __init__(self, win, size, pos, name, vertpos):
        offset = target_offsets[pos]
        super(Target, self).__init__(win,
                                     width=size, height=size,
                                     pos=[-0.36 + 0.18 * pos, vertpos + offset],
                                     lineWidth=3,
                                     name=name)

class TargetArray(object):
    def __init__(self, win, color='black', vertpos=0, drawLabels=False):
        """
        win: Parent window object
        color: Highlight color
        """
        self.targets = []
        self.win = win
        self.color = color
        for i in range(5):
            name = 'Target %d' % i
            self.targets.append(Target(win, 0.1, i, name, vertpos))
        for t in self.targets:
            t.setLineColor(COLORS['black'], colorSpace='rgb255')
        self.labels = []
        if drawLabels:
            for i, label in enumerate(['Spacebar', 'J', 'K', 'L', ';']):
                stim = visual.TextStim(self.win, text=label, pos=[-0.36 + 0.18 * i, vertpos - 0.1 + target_offsets[i]], wrapWidth=None, alignHoriz="center", height=0.045)
                self.labels.append(stim)

    def update(self, state):
        """
        state: Boolean array of target states
        """
        for i, t in enumerate(self.targets):
            if state[i]:
                #t.setLineColor(self.color, colorSpace='rgb255')
                t.setFillColor(self.color, colorSpace='rgb255')
            else:
                #t.setLineColor(COLORS['black'], colorSpace='rgb255') # black border, no fill
                #t.setFillColor(None, colorSpace='rgb255')
                t.setFillColor(BG_COLOR, colorSpace='rgb255')

    def redraw(self):
        for t in self.targets:
            t.draw()
        for t in self.labels:
            t.draw()

    def hide(self):
        for t in self.targets:
            t.setLineColor(COLORS['black'], colorSpace='rgb255')
            t.setFillColor(None, colorSpace='rgb255')
        self.win.update()

class Cue(object):
    def __init__(self, win, color, shape):
        self.win = win
        if shape == 'circle':
            self.stim = circle.Circle(win=win, fillColor=color, fillColorSpace='rgb255', lineColor=color, lineColorSpace='rgb255', radius=0.1)
        elif shape == 'triangle':
            self.stim = polygon.Polygon(win=win, edges=3, fillColor=color, fillColorSpace='rgb255', lineColor=color, lineColorSpace='rgb255', radius=0.1)
        elif shape == 'square':
            self.stim = polygon.Polygon(win=win, edges=4, fillColor=color, fillColorSpace='rgb255', lineColor=color, lineColorSpace='rgb255', radius=0.1)
        elif shape == 'pentagon':
            self.stim = polygon.Polygon(win=win, edges=5, fillColor=color, fillColorSpace='rgb255', lineColor=color, lineColorSpace='rgb255', radius=0.1)
        else:
            raise ValueError('Invalid Cue Type')

    def draw(self):
        self.stim.draw()

class Block(object):
    """
    One trial block, corresponding to a single graph
    """

    def __init__(self, experiment, block_number, subject, feedback=False, demo=False):
        """
        targets: TargetArray object
        sequence: list of states
        color: display color
        """
        self.experiment = experiment
        self.subject = subject
        self.block_number = block_number
        self.block = subject.get_block(experiment['day'], block_number)

        self.win = self.experiment['win']
        self.keyboard = self.experiment['keyboard']
        self.trial_time = self.experiment['trial_time']
        self.interkey_time = self.experiment['interkey_time']
        
        # walk: [0, 1, 2...]
        # onsets: [1,3,1,5...]
        # shape: 'circle'
        # color: RGB tuple (255)

        self.sequence = self.block['walk']

        self.trial_onsets = np.concatenate([[0], np.cumsum(self.block['SOAs'])]).astype(np.float)
        self.stimulus_ends = self.trial_onsets + self.trial_time
        self.trial_ends = np.concatenate([np.cumsum(self.block['SOAs']), 
                                          [np.sum(self.block['SOAs']) + self.trial_time]])
        self.finger_mapping = subject.finger_mapping
        self.data = []
        
        self.exp_timer = experiment['timer']
        self.block_timer = core.Clock()

        self.accuracy_record = []

        self.feedback = feedback
        self.demo = demo

    def run(self):
        io = self.experiment['io']
        keyboard = self.experiment['keyboard']
        win = self.experiment['win']

        if self.feedback:
            self.errorText = visual.TextStim(
                self.experiment['win'], text='Incorrect', font='', pos=(0.0, -0.15), height=0.045)
        self.cue = Cue(win=self.experiment['win'], color=self.block['color'], shape=self.block['shape'])
        self.targets = TargetArray(self.experiment['win'])
        self.targets.color = self.block['color']

        logging.data('block:')
        logging.data(self.block)

        self.cue.draw()
        win.update()
        io.clearEvents('all')
        keyboard.waitForPresses(keys=keypress_map.keys(), maxWait=3)

        self.trial = 0
        self.block_timer.reset()
        while self.trial < len(self.sequence):
        # while self.trial < 5:
            self.node = self.sequence[self.trial]
            target = self.finger_mapping[self.node]
            self.keytarget = target
            io.clearEvents('all')
            self.targets.update(target)
            self.draw()
            error = self.waitForResponse()
            self.showFeedback(target, error)
            self.trial += 1
        logging.flush()
        return self.accuracy_record

    def draw(self):
        #if self.trialError:
        #    self.errorText.draw()
        self.targets.redraw()
        self.win.update()

    def trialTimeLeft(self):
        """
        Each trial ends at 
        """
        return self.trial_ends[self.trial] - self.block_timer.getTime()

    def responseTimeLeft(self):
        return self.stimulus_ends[self.trial] - self.block_timer.getTime()

    def showFeedback(self, target, error):
        if self.demo:
            waitTime = 0.2
        else:
            waitTime = self.responseTimeLeft()
        if waitTime > 0:
            if self.feedback:
                if error:
                    self.targets.color = COLORS['red']
                else:
                    self.targets.color = COLORS['green']
                self.targets.update(target)
                self.draw()
                self.targets.color = self.block['color']
            core.wait(waitTime)
        if self.demo:
            waitTime = 1
        else:
            waitTime = self.trialTimeLeft()
        if waitTime > 0:
            self.win.update()
            core.wait(waitTime)

    def waitForResponse(self):
        trial = self.trial
        node = self.node
        target = self.keytarget

        trialCompleted = False
        trialError = False
        interkeyWaitTime = None
        nPrevCorrect = 0

        while not trialCompleted:
            # Set the time to wait for a keypress
            # At most the trial time remaining, if we've pressed one key then
            # up to that amount
            maxWait = self.responseTimeLeft()
            if interkeyWaitTime:
                maxWait = min(interkeyWaitTime, maxWait)

            kbe = self.keyboard.waitForPresses(keys=keypress_map.keys(), maxWait=maxWait)
            keystate = self.keyboard.state.keys()
            exp_time = self.exp_timer.getTime()
            block_time = self.block_timer.getTime()
            timeLeft = self.responseTimeLeft()

            response = [int(x in keystate) for x in response_keys]
            nCorrect = sum([x and y for x, y in zip(response, target)])
            nIncorrect = sum(response) - nCorrect
            responseIsCorrect = (response == target)
            rt = self.trial_time - timeLeft

            if responseIsCorrect:
                # Correct
                trialCompleted = True
                trialError = False
            elif (nIncorrect == 0) and (nCorrect > nPrevCorrect):
                # Pressed a new correct key
                nPrevCorrect = nCorrect
                if interkeyWaitTime is None:
                    interkeyWaitTime = self.interkey_time
            elif timeLeft < 0 and not self.demo:  # timed out on the trial
                trialCompleted = True
                trialError = True
            elif interkeyWaitTime and not kbe and not self.demo:  # interkey timeout
                trialCompleted = True
                trialError = True
            elif (nIncorrect > 0):  # incorrect response
                # wrong key
                trialCompleted = True
                trialError = True

            event_data = {'nCorrect': nCorrect,
                          'nIncorrect': nIncorrect,
                          'timedOut': timeLeft < 0,
                          'node': node,
                          'target': target,
                          'response': response,
                        #   'SOA': self.jitter[self.trial],
                          'trial': trial,
                          'blockTime': block_time,
                          'expTime': exp_time,
                          'trialTime': rt,
                          'error': trialError,
                          'completed': trialCompleted}
            self.data.append(event_data)
            logging.data(event_data)
        self.accuracy_record.append(not trialError)
        return trialError

    def save(self):
        filename = 'data/{subject}/{subject}_{day}_{block}.csv'.format(
            subject=self.subject.subject_id,
            day=self.experiment['day'],
            block=self.block_number)

        df = pd.DataFrame(self.data)
        df['subject_id'] = self.subject.subject_id
        df['block_number'] = self.block_number
        return df.to_csv(filename)

class Subject:
    """
    Graphs (networkx)
    self.graphs.lattice
    self.graphs.random
    self.graphs.modular

    Walks (node lists, 0-14)
    self.walks['modular'][1-4]
    self.walks['lattice'][1-4]
    self.walks['random'][1-4]

    Jitter
    self.SOAs['lattice'][1-4]

    Block Order
    self.blocks[day][number] = (graph, walk #) = ('modular', 1)

    Finger Mapping
    self.finger_mapping
    """
    def __init__(self, subject_id, block_order, visual_condition, n_trials):
        """
        block_order: 0 = MLML/LMLM
                     1 = LMLM/MLML

        visual_condition: 0-3 (cue=x%2, color=x//2)

        n_trials: number of trials per block
        """
        self.subject_id = subject_id
        self.block_order = block_order
        self.n_trials = n_trials
        self.shapes = cue_perms[visual_condition % 2]
        self.color = color_perms[visual_condition // 2]

        self.finger_mapping = FINGER_COMBINATIONS
        random.shuffle(self.finger_mapping)

        ## Assign the underlying graphs
        self.graphs = dict(
            demo=graphs.lattice,
            random=graphs.representative_graphs(
                graphs.generate_random_graphs(100))[0],
            modular=graphs.modular,
            lattice=graphs.lattice)

        ## Assign the graph ordering:
        day_1 = []
        day_2 = []
        day_1.append(('demo', 0))
        day_1.append(('random', 0))
        day_1.append(('random', 1))
        if self.block_order == 1:
            day_1.append(('lattice', 0))
            day_1.append(('modular', 0))
            day_1.append(('lattice', 1))
            day_1.append(('modular', 1))

            day_2.append(('modular', 2))
            day_2.append(('lattice', 2))
            day_2.append(('modular', 3))
            day_2.append(('lattice', 3))
        else:
            day_1.append(('modular', 0))
            day_1.append(('lattice', 0))
            day_1.append(('modular', 1))
            day_1.append(('lattice', 1))

            day_2.append(('lattice', 2))
            day_2.append(('modular', 2))
            day_2.append(('lattice', 3))
            day_2.append(('modular', 3))
        day_2.append(('random', 2))
        day_2.append(('random', 3))
        self.blocks = [day_1, day_2]

        self.walks = {'random': [], 'modular': [], 'lattice': [], 'demo': []}
        self.SOAs = {'random': [], 'modular': [], 'lattice': [], 'demo': []}
        for i in range(4):
            w = get_good_randomlattice_walk(self.graphs['random'], self.n_trials)
            self.walks['random'].append(w)
            self.SOAs['random'].append(get_randomlattice_SOAs(w))

            w = get_good_randomlattice_walk(self.graphs['lattice'], self.n_trials)
            self.walks['lattice'].append(w)
            self.SOAs['lattice'].append(get_randomlattice_SOAs(w))

            w = get_good_modular_walk(self.graphs['modular'], self.n_trials)
            self.walks['modular'].append(w)
            self.SOAs['modular'].append(get_modular_SOAs(w))

        w = walks.random_walk(self.graphs['demo'], 10)
        self.walks['demo'].append(w)
        self.SOAs['demo'].append(get_randomlattice_SOAs(w))

    def get_block(self, day, block):
        """
        block = {
            walk: [0,1,2...]
            SOAs: [1,3,1,5...]
            shape: 'circle'
            color: RGB tuple (255)
        }
        """
        graph, num = self.blocks[day][block]
        blockdata = {
            'walk': self.walks[graph][num],
            'SOAs': self.SOAs[graph][num],
            'shape': self.shapes[graph],
            'color': self.color[graph]
        }
        return blockdata

"""
Walk/jitter criteria:
For each walk on a modular graph, we want 
"""

"""
Walk selection code
"""

def get_clusters(walk):
    """
    Get modular graph transitions for a walk
    """
    clusters = np.zeros_like(walk)
    clusters[walk > 4] = 1
    clusters[walk > 9] = 2
    return clusters


def get_transitions(walk):
    """
    Get a boolean array of transitions
    """
    clusters = get_clusters(walk)
    transitions = np.append([0], np.diff(clusters))
    transitions = (transitions != 0)
    return transitions


def count_transitions(walk):
    """
    Get the number of transitions in a walk
    """
    return np.sum(get_transitions(walk))


def count_switchbacks(transitions):
    """
    Takes a list of trials
    0 = no transition
    1 = transition
    
    Assumes the modular graph, where two transitions in a row must
    be between the same two clusters
    """
    return np.sum((transitions[:-1]) & (transitions[1:]))


def enough_node_vists(walk):
    """
    Check that a walk visits every node at least three times
    A 100-trial block can plausibly visit every node six times
    """
    c = Counter(walk)
    for i in range(15):
        if not i in c:
            c[i] = 0

    return min(list(c.values())) > 2

def enough_transitions(walk):
    """
    Check that a walk transitions between clusters at least 9 times
    Does not include switchbacks
    """
    n_transitions = count_transitions(walk)
    t = get_transitions(walk)
    n_switchbacks = count_switchbacks(t)
    return n_transitions - n_switchbacks > 8

def get_good_modular_walk(graph, n):
    """
    Finds a modular graph walk that:
    visits every node at least 3 times
    transitions between clsuters at least 9 times
    """
    walk = np.array(walks.random_walk(graph, n))
    while not (enough_node_vists(walk) and enough_transitions(walk)):
        walk = np.array(walks.random_walk(graph, n))
    return walk

def get_good_randomlattice_walk(graph, n):
    """
    Finds a walk that:
    visits every node at least 3 times
    """
    walk = np.array(walks.random_walk(graph, n))
    while not enough_node_vists(walk):
        walk = np.array(walks.random_walk(graph, n))
    return walk

"""
Efficiency code
"""

from scipy.io import loadmat
# Using the HRF from SPM
hrf_25 = loadmat('hrf_25.mat')
hrf_25 = hrf_25['hrf_25'].flatten()

def calculate_modular_efficiency(transitions, SOAs, effect_size):
    """
    transitions: list of where a walk transitions betwen clusters
    SOAs: List of stimulus onset times
    Expected effect size difference
    """
    t_within = np.zeros(sum(SOAs)+1)
    t_between = np.zeros(sum(SOAs)+1)
    t_within[np.cumsum(SOAs)[~transitions[:-1]]] = 1
    t_between[np.cumsum(SOAs)[transitions[:-1]]] = 1

    t_between *= effect_size

    # Upsample, convolve with the HRF, and downsample
    r_within = np.repeat(t_within, 4)
    r_within = np.convolve(r_within, hrf_25)
    r_within = r_within[::4]

    r_between = np.repeat(t_between, 4)
    r_between = np.convolve(r_between, hrf_25)
    r_between = r_between[::4]

    r_total = r_within + r_between

    # D = [X; X2]';
    # c = [0 1];
    # 1/(c*pinv(D'*D)*c')
    # corr(X', X2')

    # eff = 1/(c@pinv(D@D.T)@c.T
    D = np.vstack((r_within, r_between))
    c = np.array([0, 1])
    #eff = 1/(c@pinv(D@D.T)@c.T)
    eff = 1 / matmul(matmul(c, pinv(matmul(D, D.T))), c.T)
    
    return {
        'transitions': transitions,
        't_between': t_between,
        't_within': t_within,
        'r_between': r_between,
        'r_within': r_within,
        'eff': eff
    }


def get_modular_SOAs(walk):
    n_trials = len(walk)
    clusters = get_clusters(walk)
    transitions = get_transitions(walk)

    efficiency_list = []
    for i in range(200):
        SOAs = np.repeat([2,4,6],int(np.ceil(n_trials/3)))
        SOAs = SOAs[:(n_trials-1)]
        random.shuffle(SOAs)
        SOAs = np.array(SOAs)

        t = calculate_modular_efficiency(transitions, SOAs, 1.1)
        t['SOAs'] = SOAs
        efficiency_list.append(t)
    
    df = pd.DataFrame(efficiency_list)
    df = df.sort_values('eff', ascending=False).reset_index(drop=True)

    return df.iloc[0].SOAs

def get_randomlattice_SOAs(walk):
    n_trials = len(walk)
    SOAs = np.repeat([2, 4, 6], int(np.ceil(n_trials / 3)))
    SOAs = SOAs[:(n_trials-1)]
    random.shuffle(SOAs)
    SOAs = np.array(SOAs)
    return SOAs


"""
Tests
"""
def test_modular_walk():
    print("Testing modular walk generation")
    graph = graphs.modular
    for i in range(10):
        walk = get_good_modular_walk(graph, 100)
        n_transitions = count_transitions(walk)
        t = get_transitions(walk)
        n_switchbacks = count_switchbacks(t)
    
        c = Counter(walk)
        for j in range(15):
            if not j in c:
                c[j] = 0

        min_visits = min(list(c.values()))
        
        print("Walk %d: %d transitions and %d minimum visits" % (i, n_transitions - n_switchbacks, min_visits))
        assert n_transitions - n_switchbacks > 8
        assert min_visits > 2

    print("Generated 10 good walks")
    print("Modular walk generation OK\n")


def test_modular_SOAs():
    print("Testing modular SOA generation")
    graph = graphs.modular
    walk = get_good_modular_walk(graph, 100)
    for i in range(10):
        SOAs = get_modular_SOAs(walk)
        transitions = get_transitions(walk)
        eff = calculate_modular_efficiency(transitions, SOAs, 1.1)
        print("efficiency: %f" % eff['eff'])
    print(SOAs)
    print("Modular SOA generation OK\n")

def test_random_walk():
    print("Testing random walk generation")
    graph = graphs.representative_graphs(
        graphs.generate_random_graphs(100))[0]
    for i in range(10):
        walk = get_good_randomlattice_walk(graph, 100)
    
        c = Counter(walk)
        for j in range(15):
            if not j in c:
                c[j] = 0

        min_visits = min(list(c.values()))
        print(walk[:10])
        print("Walk %d: %d minimum visits" % (i, min_visits))
        assert min_visits > 2

    print("Generated 10 good walks")
    print("Random walk generation OK\n")

def test_random_SOAs():
    print("Testing random SOA generation...")
    graph = graphs.representative_graphs(
        graphs.generate_random_graphs(100))[0]
    walk = get_good_randomlattice_walk(graph, 100)
    SOAs = get_randomlattice_SOAs(walk)
    print("Random SOA generation OK\n")

def test_subject_generation():
    print("Testing subject initialization")
    s = Subject('test', 0, 3, 100)
    s = Subject('test', 1, 3, 100)
    print(s.blocks)
    print(s.color)
    print(s.shapes)
    print("Subject initialization ok\n")


def test_subject_blocks():
    print("Testing subject block retrieval")
    s = Subject('test', 0, 3, 100)
    for day in range(2):
        for block in range(6):
            b = s.get_block(day,block)
    print(b)
    print("Subject block retrieval ok\n")

def test_block_creation():
    print("Testing block creation")
    s = Subject('test', 0, 3, 100)
    experiment = {
        'win': None,
        'io': None,
        'keyboard': None,
        'day': 0,
        'timer': None,
        'trial_time': 1,
        'interkey_time': 1
    }
    b = Block(experiment,0,s)
    assert len(b.trial_onsets) == len(b.sequence)
    assert len(b.trial_ends) == len(b.sequence)
    assert len(b.stimulus_ends) == len(b.sequence)
    assert b.stimulus_ends[-1] == b.trial_ends[-1]
    assert b.trial_onsets[0] == 0
    assert b.stimulus_ends[0] == b.trial_time
    print("Demo times:")
    print(" Trial Onsets:", b.trial_onsets)
    print("Stimulus Ends:", b.stimulus_ends)
    print("   Trial Ends:", b.trial_ends.astype(np.float))
    print("Block creation OK\n") 


# def test_demo():
#     print("Testing a short demo block")
#     s = Subject('test', 0, 3, 20)
#     win = visual.Window((1920,1000), monitor="testMonitor", units='height', fullscr=False)
#     win.setColor(BG_COLOR, colorSpace='rgb255')
#     expTimer = core.Clock()
#     win.flip()

#     block = Block(experiment=experiment,
#                     block_number=0,
#                     subject=subject,
#                     feedback=False)
#     block.run()
#     # block.save()

#     win.close()

def run_all_tests():
    test_modular_walk()
    test_modular_SOAs()

    test_random_walk()
    test_random_SOAs()

    test_subject_generation()

    test_block_creation()

    # test_block_run()
