# March-Madness-2020
Our prediction model for March Madness 2020


[**Full Presentation Here**](https://github.com/ngocdinh1410/March-Madness-2020/blob/master/GUCCIGANG%20Report.pptx)


[**Full R file here**](https://github.com/ngocdinh1410/March-Madness-2020/blob/master/final.R)

[**Our tableau visualization**](https://github.com/ngocdinh1410/March-Madness-2020/blob/master/March%20Madness%20Test%20Data%20Visualization.twbx)


**<h2>Process to tackle the modeling</h2>**


![alt text](https://github.com/ngocdinh1410/March-Madness-2020/blob/master/Process.PNG "Our process")
<h2>Data Pre-processing:</h2>
<h3>Derived Novel Variables</h3>
<li>Team’s winning rate – Calculated as ‘wins/(wins+losses)’.</li>
<li>Team’s win-loss ratio – Calculated as ‘wins/losses’.</li>
<li>Coach’s winning rate – Calculated as ‘wins/(wins+losses)’.</li>
<li>Coach’s win-loss ratio – Calculated as ‘wins/losses’.</li>
<li>Seed Differences – Calculated as ‘strong seed - weak seed’.</li>
<li>Whether team 1 wins – Derived by scores after randomly switch team 1 and team 2.</li>
<h3>Variables Removal</h3>
<li>Remove two teams’s variables in association with wins and losses.</li>
<h3>Split Data (Used in prediction model building)</h3>
<li>Training Data: 2002 – 2018; Test Data: 2019</li>
<li>Data needed prediction: 2020</li>

<h2> Visualization</h2>

<h3>Assist Scoring:</h3>


![alt text](https://github.com/ngocdinh1410/March-Madness-2020/blob/master/Assist%20scoring%20March%20Madness.png "picture 1")


<h3>Dashboard of team 1 vs team 2 stats:</h3>


![alt text](https://github.com/ngocdinh1410/March-Madness-2020/blob/master/Dashboard%201%20MARCH%20MADNESS.png "picture 2")


<h3>Does seed matters?:</h3>


![alt text](https://github.com/ngocdinh1410/March-Madness-2020/blob/master/coach%20vs%20season%20march%20madness.png "picture 3")


<h2>Modeling:</h2>

<li>Logistic regression</li>
<li>Random forest</li>
<li>Linear discriminant analysis (LDA)</li>
<li>Quadratic discriminant analysis (QDA)</li>
<li>Support vector machine (SVM)</li>

<h2>Model Result:</h2>


![alt text](https://github.com/ngocdinh1410/March-Madness-2020/blob/master/Model%20result.PNG "Model Result")
