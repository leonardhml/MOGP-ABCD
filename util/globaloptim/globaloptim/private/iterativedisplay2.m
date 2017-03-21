function textArea = iterativedisplay2(text, title, x, y, w, h) 
%ITERATIVEDISPLAY2 displays text in a scrollable text area
%
%  TEXTAREA = ITERATIVEDISPLAY2 returns a handle to the text area which
%  can use Java TextArea methods such as setText or append.
%
%  TEXTAREA = ITERATIVEDISPLAY2(TEXT) displays TEXT at startup.
%
%  TEXTAREA = ITERATIVEDISPLAY2(TEXT, TITLE) sets the title of 
%  the window to TITLE.
%
%  TEXTAREA = ITERATIVEDISPLAY2(TEXT, TITLE, x, y, w, h) sets the
%  location and size of the window.

%   Copyright 2003-2010 The MathWorks, Inc.

import com.mathworks.mwswing.*;
import javax.swing.*;
import javax.swing.border.*
import java.awt.*

if nargin < 6
    h = 300;
end
if nargin < 5
    w = 650;
end
if nargin < 4
    y = 400;
end  
if nargin < 3  
    x = 400;
end
if nargin < 2  
    title = '';
end
if nargin < 1  
   text = ''; 
end
    
idFrame = MJFrame(title);
idFrame.setBounds(x, y, w, h);
cp = idFrame.getContentPane;
cp.setLayout(BorderLayout);
textArea = MJTextArea;
textArea.setEditable(false);
textArea.setText(text);
taSP = MJScrollPane(textArea);
cp.add(taSP);
idFrameHandle = handle(idFrame,'callbackProperties'); 
set(idFrameHandle,'WindowClosingCallback', {@closeWindow, idFrame, textArea});
idFrame.setVisible(true);

%--------------------------------------------------------------------------
function closeWindow(obj, evd, h, ta)
if nargin > 2 && ~isempty(h) && ishandle(h)
    h.dispose;
end
