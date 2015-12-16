#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import *
from tkinter.ttk import *
from tkinter.font import Font

import simpcalc


def autoclose(expr):
    count = 0
    for ch in expr:
        if ch == '(':
            count += 1
            continue
        elif ch == ')':
            if count:
                count -= 1
    return expr + ')' * count


class CopiableText(Text):

    def __init__(self, master, **kw):
        Text.__init__(self, master, **kw)
        self.bind('<Control-c>', self.copy)
        self.bind('<Control-x>', self.cut)

    def copy(self, event=None):
        self.clipboard_clear()
        text = self.get("sel.first", "sel.last")
        self.clipboard_append(text)

    def cut(self, event):
        self.copy()
        self.delete("sel.first", "sel.last")


class App(Frame):

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.calculator = simpcalc.Calculator('ans', True)
        self.returned = False

        self.buttongrid = (
            ('7', '8', '9', '÷', 'C', '⌫'),
            ('4', '5', '6', '×', '(', ')'),
            ('1', '2', '3', '-', '^', '√'),
            ('0', '.', 'e', '+', '=')
        )
        self.master.geometry('320x320')
        self.master.title('Calculator')
        self.master.wm_resizable(width=False, height=False)
        self.createWidgets()

    def createWidgets(self):
        top = self.winfo_toplevel()
        font = Font(size=12)
        fontb = Font(size=12, weight='bold')
        self.grid(sticky=NSEW)
        top.grid_rowconfigure(0, weight=1)
        top.grid_columnconfigure(0, weight=1)
        for i in range(6):
            self.columnconfigure(i, minsize=40, weight=1)
        self.rowconfigure(0, minsize=50, weight=4)
        self.rowconfigure(1, minsize=25, weight=2)
        for i in range(2, 6):
            self.rowconfigure(i, minsize=25, weight=1)

        self.text_history = CopiableText(self, height=5, font=font)
        self.text_history.grid(row=0, column=0, columnspan=6, sticky=NSEW)
        self.text_history.tag_configure("bold", font=fontb)
        self.text_input = CopiableText(self, height=2, font=font)
        self.text_input.grid(row=1, column=0, columnspan=6, sticky=NSEW)
        self.text_input.tag_configure("right", justify='right')
        self.text_input.tag_configure("bold", font=fontb, justify='right')
        self.text_input.tag_configure("err", foreground='red', justify='right')

        self.buttons = {}
        for y, row in enumerate(self.buttongrid):
            for x, t in enumerate(row):
                self.buttons[t] = Button(
                    self, text=t, command=self.input_text(t))
                self.buttons[t].grid(
                    row=y + 2, column=x, sticky=NSEW, padx=2, pady=2)

        self.buttons['='].grid(columnspan=2)
        self.buttons['=']['default'] = 'active'
        self.buttons['=']['command'] = self.submit
        self.buttons['⌫']['command'] = self.backspace
        self.buttons['C']['command'] = self.clear

        keys = '!%()+,-.0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ\\^abcdefghijklmnopqrstuvwxyz'
        for i in keys:
            self.bind_all(i, self.input_text(i))
        for i in '0123456789':
            self.bind_all('<KP_%s>' % i, self.input_text(i))
        self.bind_all('<KP_Decimal>', self.input_text('.'))
        self.bind_all('<KP_Add>', self.input_text('+'))
        self.bind_all('<KP_Subtract>', self.input_text('-'))
        self.bind_all('<KP_Multiply>', self.input_text('×'))
        self.bind_all('<KP_Divide>', self.input_text('÷'))
        self.bind_all('*', self.input_text('×'))
        self.bind_all('/', self.input_text('÷'))
        self.bind_all('=', self.submit)
        self.bind_all('<Delete>', self.clear)
        self.bind_all('<Return>', self.submit)
        self.bind_all('<KP_Enter>', self.submit)
        self.bind_all('<BackSpace>', self.backspace)
        self.bind_all('<Control-v>', self.paste)

    def input_text(self, t):
        def input_char(event=None):
            if not (event and event.widget is self.text_input):
                if self.returned and t not in '!%+-.\\^×÷':
                    self.clear()
                self.returned = False
                self.text_input.insert(INSERT, t)
                self.text_input.tag_add("right", '1.0', END)
        return input_char

    def backspace(self, event=None):
        if not (event and event.widget is self.text_input):
            self.text_input.delete('insert -1 chars')

    def clear(self, event=None):
        if not (event and event.widget is self.text_input):
            self.text_input.delete('1.0', END)

    def paste(self, event=None):
        if not (event and event.widget is self.text_input):
            text = self.text_input.selection_get(selection='CLIPBOARD')
            self.text_input.insert('insert', text)

    def submit(self, event=None):
        expr = self.text_input.get('1.0', END)
        expr = expr.strip().replace('\r', '').replace('\n', '')
        if not expr:
            return
        try:
            r = self.calculator.eval(expr)
        except simpcalc.CalculatorError as ex:
            self.text_input.tag_add("err", '1.%d' %
                                    ex.pos, '1.%d' % (ex.pos + ex.length))
            return
        ret = self.calculator.format(r)
        self.text_input.replace('1.0', END, ret)
        self.text_input.tag_add("bold", '1.0', 'end -1 chars')
        self.text_input.mark_set(INSERT, END)
        self.text_history.insert(END, '%s = %s\n' % (autoclose(expr), ret))
        self.text_history.tag_add("bold", 'end -%d chars' %
                                  (len(ret) + 5), 'end -2 chars')
        self.text_history.mark_set(INSERT, END)
        self.text_history.see(END)
        self.returned = True

style = Style()
style.theme_use('clam')

app = App()
app.mainloop()
