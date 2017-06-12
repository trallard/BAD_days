---
layout: page
permalink: /1-outline/
title: Day 1 outline
category: module
---


**Description**

The Bioinformatics Awareness Days are days devoted to Bioinformatics. As a Marie S. Curie postdoctoral fellow @ SITRaN and Computer Science Departments at the University of Shegffield, I organized such an event to give my practical contribution in Bioinformatics, although I have to state that all this started in Naples at the Telethon Institute of Genetics and Medicine, with the great Bionformatics Core.
Together with Tania Allard and Micke Croucer at the DCS, we decided to publicly divulgate this material to all the interested scientific community.
The sessionsare self contained and a full run should last at most 2 hours.

Prerequisites: Addressed to researchers already familiar with R  and basics of statistics

Day1 focus: "R&Bioconductor and Multiple hypothesis testing"

Program:


- R&Bioconductor
- Basic data types, Classes and Methods
- Useful R packages and functions
- Multiple hypothesis Testing (FDR, Bonferroni and Benjamini correction)
- Examples (exercise notebook)

 




## Outline
This tutorial session will provide a brief introduction to Data Science and data manipulation using R and the Jupyter notebooks.

The hands-on components of this module are contained in the following tutorials:

<ul >
{% assign pages_list = site.pages %}
{% for node in pages_list %}
{% if node.title != null %}
{% if node.layout == "default" %}
{%if node.tags %}

{% for tag in node.tags %}
{% if tag == 'Day1' %}
<!-- Note you need to prepend the site.baseurl always-->
<li><a href="{{site.baseurl}}{{ node.url }}">{{ node.title }}</a>
</li>
{% endif %}
{% endfor %}

{% endif %}
{%endif%}
{% endif %}
{% endfor %}
</ul>
