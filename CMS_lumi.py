# uncompyle6 version 3.9.0
# Python bytecode version base 2.7 (62211)
# Decompiled from: Python 3.7.1 (default, Dec 14 2018, 13:28:58) 
# [Clang 4.0.1 (tags/RELEASE_401/final)]
# Embedded file name: /work/jinw/CMSSW_10_2_15_patch2/src/cleanup_2Dunfold/CMS_lumi.py
# Compiled at: 2023-03-06 18:58:44
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
cmsText = 'CMS'
cmsTextFont = 61
writeExtraText = True
extraText = 'Preliminary'
extraTextFont = 52
lumiTextSize = 0.65
lumiTextOffset = 0.03
cmsTextSize = 0.95
cmsTextOffset = 0.5
relPosX = 0.645
relPosX = 0.145
relPosY = -0.011
relExtraDY = 1.2
extraOverCmsTextSize = 0.76
lumi_13TeV = '20.1 fb^{-1}'
lumi_8TeV = '19.7 fb^{-1}'
lumi_7TeV = '5.1 fb^{-1}'
lumi_sqrtS = ''
drawLogo = False

def CMS_lumi(pad, iPeriod, iPosX):
    outOfFrame = False
    if iPosX / 10 == 0:
        outOfFrame = True
    if iPosX / 10 == 0:
        alignX_ = 1
    if iPosX == 0:
        alignY_ = 1
    if iPosX / 10 == 1:
        alignX_ = 1
    if iPosX / 10 == 2:
        alignX_ = 2
    if iPosX / 10 == 3:
        alignX_ = 3
    alignY_ = 1
    alignX_ = 1
    align_ = 10 * alignX_ + alignY_
    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025
    pad.cd()
    lumiText = ''
    if iPeriod == 1:
        lumiText += lumi_7TeV
        lumiText += ' (7 TeV)'
    elif iPeriod == 2:
        lumiText += lumi_8TeV
        lumiText += ' (8 TeV)'
    elif iPeriod == 3:
        lumiText = lumi_8TeV
        lumiText += ' (8 TeV)'
        lumiText += ' + '
        lumiText += lumi_7TeV
        lumiText += ' (7 TeV)'
    elif iPeriod == 4:
        lumiText += lumi_13TeV
        lumiText += ' (13 TeV)'
    elif iPeriod == 7:
        if outOfFrame:
            lumiText += '#scale[0.85]{'
        lumiText += lumi_13TeV
        lumiText += ' (13 TeV)'
        lumiText += ' + '
        lumiText += lumi_8TeV
        lumiText += ' (8 TeV)'
        lumiText += ' + '
        lumiText += lumi_7TeV
        lumiText += ' (7 TeV)'
        if outOfFrame:
            lumiText += '}'
    elif iPeriod == 12:
        lumiText += '8 TeV'
    elif iPeriod == 0:
        lumiText += lumi_sqrtS
    print lumiText
    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)
    extraTextSize = extraOverCmsTextSize * cmsTextSize
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(lumiTextSize * t)
    latex.DrawLatex(1 - r, 1 - t + lumiTextOffset * t, lumiText)
    if outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(cmsTextSize * t)
        latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText)
    pad.cd()
    posX_ = 0
    if iPosX % 10 <= 1:
        posX_ = l + relPosX * (1 - l - r)
    elif iPosX % 10 == 2:
        posX_ = l + 0.5 * (1 - l - r)
    elif iPosX % 10 == 3:
        posX_ = 1 - r - relPosX * (1 - l - r)
    posY_ = 1 - t + lumiTextOffset * t
    if not outOfFrame:
        if drawLogo:
            posX_ = l + 0.045 * (1 - l - r) * W / H
            posY_ = 1 - t - 0.045 * (1 - t - b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15 * H / W
            yl_1 = posY_
            CMS_logo = rt.TASImage('CMS-BW-label.png')
            pad_logo = rt.TPad('logo', 'logo', xl_0, yl_0, xl_1, yl_1)
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw('X')
            pad_logo.Modified()
            pad.cd()
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize * t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if writeExtraText:
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize * t)
                latex.DrawLatex(posX_, posY_ - relExtraDY * cmsTextSize * t, extraText)
    elif writeExtraText:
        if iPosX == 0:
            posX_ = l + relPosX * (1 - l - r)
            posY_ = 1 - t + lumiTextOffset * t
        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize * t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)
    pad.Update()