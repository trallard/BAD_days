---
layout: page
permalink: /2-outline/
title: Day 2 outline
category: module
---


**Description**
This second session is self contained and an average expert user would take about two hours for a full run 

Prerequisites: Addressed to researchers already familiar with R and basics of statistics

Day 2 Program:
- Dimensionality Reduction Methods
- Analysis of Variance and Generalized Linear Models
- Case studies and Examples in R
---

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
{% if tag == 'Day2' %}
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
