# Psychopy on windows seems to want to be imported first
from psychopy import visual, core, event, logging, gui
from psychopy.visual import rect
from psychopy.visual import polygon
from psychopy.visual import ImageStim
from psychopy.visual import TextStim
from psychopy.visual import circle
# from psychopy.visual.shape import BaseShapeStim
from psychopy.tools.attributetools import attributeSetter, setAttribute

from psychopy.iohub import launchHubServer

# just want to make sure this exists
# Otherwise iohub gives an ambiguous warning
import msgpack_numpy

import numpy as np
import pandas as pd
import csv
import sys
from datetime import datetime
import os
import errno
import random
import pickle
from networkx.readwrite import json_graph
import json
import shutil

import sys
sys.path.append('.')
import graphs
import walks

from models import Subject, TargetArray, Block, BG_COLOR, COLORS

io = launchHubServer()
keyboard = io.devices.keyboard

BLOCK_LENGTH = 100
#PRACTICE_LENGTH = 10

PRACTICE_TRIAL_TIME = 6
TRIAL_TIME = 1.0

BREAK_TIME = 30 # Tme between blocks

INTERKEY_TIME = 0.2 # seconds between keys
SHOW_INSTRUCTIONS = True
SHOW_QUIZ = False
FULLSCREEN = False
RESOLUTION = [1920, 1200]

class Instructions(object):
    def __init__(self, win):
        self.win = win

    def makeTextStim(self, text, vertpos=0):
        return TextStim(
            self.win,
            text=text,
            wrapWidth=1.0,
            pos=(0.0, vertpos),
            alignHoriz="center",
            height=0.08,
            color="#FFFFFF",
            units="norm")

    def run(self):
        pressToContinue = TextStim(self.win, text="<press any key to continue>", pos=(0, -0.9), height=0.06, units="norm", alignHoriz="center")
        #fixation = TextStim(win, text="+", height=0.1, color="#FFFFFF")

        instruct_1_text = """You will see five squares on the screen, which correspond with 'Spacebar', 'J', 'K', 'L', and ';'."""
        instruct1 = self.makeTextStim(instruct_1_text, vertpos=0.3)

        # example_targets_blank = TargetArray(win, vertpos=-0.1, drawLabels=False)
        example_targets_labeled = TargetArray(win, vertpos=-0.1, drawLabels=True, color=COLORS['purple'])
        

        instruct_2_text = "When a square lights up, your job is to press that key as quickly as possible.\n\n"
        instruct_2_text += "If two squares light up, press them at the same time (or as close as you can)."
        instruct2 = self.makeTextStim(instruct_2_text, vertpos=0.3)

        instruct_3_text = "You should press 'SPACEBAR' with your thumb, 'J' with your index finger, 'K' with your middle finger, 'L' with your ring finger, and ';' with your pinky.\n\n"
        instruct_3_text += "Hold your thumb over space and one finger on each of the subsequent letters (See image below).\n\n"
        instruct_3_text += "If you need to adjust where the keyboard is for this to be comfortable, please do so now."
        instruct3 = self.makeTextStim(instruct_3_text, vertpos=0.4)
        handStim = ImageStim(self.win, image='img/hand.png', size=(0.6, 0.45), pos=(0, -0.3))

        instruct_4_text = "You will have a limited amount of time to make your response. "
        instruct_4_text += "The squares will remain on the screen for about 1 second, followed by a blank screen for another brief period.\n\n"

        instruct_4_text += "To answer correctly, press the correct keys before the squares disappear. "
        instruct_4_text += "If you press the wrong keys, the trial will count as a mistake. "
        instruct_4_text += "If you take too long to answer, the trial will count as a mistake.\n\n"

        instruct_4_text += "This may be difficult at first, but will become easier with practice. "
        instruct_4_text += "In either case, the blank screen will advance to the next trial after a few seconds."
        instruct4 = self.makeTextStim(instruct_4_text, vertpos=0)

        instruct_5_text = "You will complete six sequences of trials. Each sequence will last around 6-7 minutes. "
        instruct_5_text += "You can take a short break after each sequence.\n\n"
        instruct_5_text += "There are three distinct types of sequences. You may or may not notice differences between them. \n\n"
        instruct_5_text += "Each type of sequence will show a different shape in the middle of the screen for a few seconds before the sequence starts. The squares will also show up as a different color for each."
        instruct5 = self.makeTextStim(instruct_5_text, vertpos=0)

        instruct_6_text = "Once you finish the experiment, your accuracy will be calculated as the percentage of trials without a mistake. "
        instruct_6_text += "If this is at least 90%, you will receive a payment bonus of $5."
        instruct6 = self.makeTextStim(instruct_6_text, vertpos=0) 
        
        example_targets_labeled.redraw()
        instruct1.draw()
        pressToContinue.draw()
        win.flip()
        event.waitKeys() 

        instruct2.draw()
        example_targets_labeled.update([0,1,0,1,0])
        example_targets_labeled.redraw()
        win.flip()
        event.waitKeys()
        
        instruct3.draw()
        handStim.draw()
        win.flip()
        event.waitKeys()

        instruct4.draw()
        win.flip()
        event.waitKeys()
        
        instruct5.draw()
        win.flip()
        event.waitKeys()

        instruct6.draw()
        win.flip()
        event.waitKeys()

    def day2(self):
        day2_text = "The task for this session is similar to yesterday. "
        day2_text += "You will respond six additional sets of trials. "

        pressToContinue = TextStim(self.win, text="<press any key to continue>", pos=(0, -0.9), height=0.06, units="norm", alignHoriz="center")

        ready = self.makeTextStim(day2_text)
        ready.draw()
        pressToContinue.draw()
        win.flip()
        event.waitKeys()

    def ready_practice(self):
        ready_practice_text = "You're almost ready to begin! "
        ready_practice_text += "To help you with the main task, you'll first have a chance to practice. "
        ready_practice_text += "This will take around a minute and will not count towards your accuracy bonus.\n\n"
        ready_practice_text += "Trials will be slower for the practice. Also, unlike the main task, you will know if you answered correctly."
        ready_practice_text += "Squares will change to green if you answered correctly, and red if incorrectly. "
        ready_practice_text += "THIS WILL NOT HAPPEN AFTER THE PRACTICE.\n\n"
        ready_practice_text += "Ready? Press any button!"

        ready = self.makeTextStim(ready_practice_text)
        ready.draw()
        win.flip()
        event.waitKeys()

    def ready_start(self):
        ready_start_text = "The first sequence will begin now. Ready? Press any button!"

        ready = self.makeTextStim(ready_start_text)
        ready.draw()
        win.flip()
        event.waitKeys()

    def ready_n(self, n):
        if n == 0:
            self.ready_start()
        else:
            ready_n_text = "Great job! You've finished %d of 6 sets of sequences. "
            ready_n_text += "You can take up to 30 seconds break if you need to. When you're ready to continue press any button."

            ready = self.makeTextStim(ready_n_text % n)
            ready.draw()
            win.flip()
            event.waitKeys(maxWait=BREAK_TIME)

    def finished(self, percentage_correct):
        ready = self.makeTextStim("Great job, you're done!\n Your accuracy was %0.1f%%" % percentage_correct)
        ready.draw()
        win.flip()
        event.waitKeys()

class Quiz(object):
    def __init__(self, win, quizFile):
        self.win = win
        self.quizQ = visual.TextStim(self.win, text="blank", pos=(0.0, 0.3), wrapWidth=1, alignHoriz="center", height=0.055)
        self.quizOpt1 = visual.TextStim(self.win, text="blank", pos=(-0.5, 0.15), wrapWidth=None, alignHoriz="left", height=0.045)
        self.quizOpt2 = visual.TextStim(self.win, text="blank", pos=(-0.5, 0.0), wrapWidth=None, alignHoriz="left", height=0.045)
        self.quizOpt3 = visual.TextStim(self.win, text="blank", pos=(-0.5, -0.15), wrapWidth=None, alignHoriz="left", height=0.045)
        self.quizResult = visual.TextStim(self.win, text="blank", pos=(0.0, 0.0), wrapWidth=1, alignHoriz="center", height=0.045)
        self.quizList = pd.read_csv(quizFile)

    def run(self, quizNumber):
        nIncorrect = 0
        trials = self.quizList[self.quizList.quiz == quizNumber]
        for trial in trials.itertuples():
            self.quizQ.setText(trial.question)
            self.quizOpt1.setText("1) " + trial.answer_1)
            self.quizOpt2.setText("2) " + trial.answer_2)
            self.quizOpt3.setText("3) " + trial.answer_3)

            self.drawQuiz()

            key = event.waitKeys(keyList=('1', '2', '3'))
            response = int(key[0][0])
            correct = (response == trial.correct_ans)
            if not correct:
                nIncorrect += 1
            logging.info('Quiz Question: %s' % trial.question)
            logging.info('Quiz Answer 1: %s' % trial.answer_1)
            logging.info('Quiz Answer 2: %s' % trial.answer_2)
            logging.info('Quiz Answer 3: %s' % trial.answer_3)
            logging.info('Quiz Correct Answer: %d' % trial.correct_ans)
            logging.info('Quiz Response: %d' % response)
            logging.flush()
        if nIncorrect > 0:
            if nIncorrect == 1:
                self.quizResult.setText('You answered 1 question wrong. Press any key to review the instructions.')
            else:
                self.quizResult.setText('You answered %d questions wrong. Press any key to review the instructions.' % nIncorrect)
            self.drawResult()
            event.waitKeys()
            return False
        else:
            self.quizResult.setText('You answered all questions correctly! Press any key to move on to the next phase.')
            self.drawResult()
            event.waitKeys()
            return True

    def drawQuiz(self):
        self.quizQ.draw()
        self.quizOpt1.draw()
        self.quizOpt2.draw()
        self.quizOpt3.draw()
        self.win.flip()

    def drawResult(self):
        self.quizResult.draw()
        self.win.flip()

def get_day():
    dayDlg = gui.Dlg(title="Image Viewing Study", pos=(-1.5, 0))
    dayDlg.addField('Day','0')
    dayDlg.show()

    if gui.OK:
        day = int(dayDlg.data[0])
    else:
        sys.exit()

    if day < -1 or day > 1:
        raise ValueError

    return day

def new_subject():
    # get subjID
    subjDlg = gui.Dlg(title="Image Viewing Study", pos=(-1.5,0))
    subjDlg.addField('Subject ID:','demo')
    subjDlg.addField('Block Order:',0)
    subjDlg.addField('Condition:', 0)
    subjDlg.show()

    if gui.OK:
        subject_id=subjDlg.data[0]
        block_order = int(subjDlg.data[1])
        visual_condition = int(subjDlg.data[2])
    else:
        sys.exit()

    subject_id ='%s_%s' % (subject_id, datetime.today().strftime('%Y%m%d_%H%M%S'))

    subject = Subject(subject_id, block_order, visual_condition, BLOCK_LENGTH)

    return subject

def load_subject():
    filename = gui.fileOpenDlg(tryFilePath='data',allowed='*.pkl')[0]
    with open(filename, 'rb') as fp:
        subject = pickle.load(fp)
    return subject

def save_subject(subject):
    subject_id = subject.subject_id
    filename = 'data/%s/%s.pkl' % (subject_id, subject_id)

    try:
        os.makedirs('data/%s' % subject_id)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    with open(filename, 'wb') as fp:
        pickle.dump(subject, fp)

if __name__ == "__main__":
    logging.console.setLevel(logging.DATA)
    win = visual.Window(RESOLUTION, monitor="testMonitor", units='height', fullscr=False)
    win.setColor(BG_COLOR, colorSpace='rgb255')
    targets = TargetArray(win)
    targetsLabeled = TargetArray(win, drawLabels=True)
    expTimer = core.Clock()
    win.flip()

    day = get_day()

    if day == -1:
        from models import run_all_tests
        run_all_tests()
        raise Exception
    
    if day == 0:
        subject = new_subject()
        save_subject(subject)
        subject_id = subject.subject_id

        # ## Copy current version of the script
        scriptname = os.path.basename(__file__)
        shutil.copyfile(scriptname, 'data/%s/experiment_code_backup.py' % subject_id)
        shutil.copyfile('models.py', 'data/%s/models_backup.py' % subject_id)

    else:
        subject = load_subject()
        subject_id = subject.subject_id

    if FULLSCREEN:
        win.winHandle.maximize()
        win.winHandle.activate()
        win.fullscr = True
        win.winHandle.set_fullscreen(True)
        win.flip()

    quiz = Quiz(win, 'quiz.csv')
    instructions = Instructions(win)

    ## Log current experiment setup
    logging.LogFile(f='data/%s/%s_day_%d.log' % (subject_id, subject_id, day), level=logging.DATA)
    logging.data('Subject ID is %s' % subject_id)
    logging.data('Day is %d' % day)
    logging.flush()

    accuracyRecord = []

    experiment = {
        'win': win,
        'io': io,
        'keyboard': keyboard,
        'day': day,
        'timer': expTimer,
        'trial_time': TRIAL_TIME,
        'interkey_time': INTERKEY_TIME
    }

    if day == 0:
        """ Initial instructions """
        if SHOW_QUIZ:
            finishedQuiz = False
            while not finishedQuiz:
                instructions.run()
                finishedQuiz = quiz.run(1)
        elif SHOW_INSTRUCTIONS:
            instructions.run()

        """ Trial Block """
        instructions.ready_practice()
        block = Block(experiment=experiment,
                      block_number=0,
                      subject=subject,
                      feedback=True)
        block.run()
        block.save()

        for i in range(1,7):
            instructions.ready_n(i-1)
            block = Block(experiment=experiment,
                        block_number=i,
                        subject=subject,
                        feedback=False)
            accuracyRecord.append(block.run())
            block.save()
    else:
        instructions.day2()
        for i in range(6):
            instructions.ready_n(i)
            block = Block(experiment=experiment,
                        block_number=i,
                        subject=subject,
                        feedback=False)
            accuracyRecord.append(block.run())
            block.save() 

    trial_scores = [np.mean(x)*100 for x in accuracyRecord]
    total_score = np.mean(trial_scores)

    logging.data('Trial accuracy scores are: %s' % trial_scores)
    logging.data('Total accuracy is: %s' % total_score)
    logging.flush()

    instructions.finished(total_score)
