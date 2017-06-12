---
layout: page
permalink: /1-outline/
title: Day 1
category: module
---

**Prerequisites:** Addressed to researchers already familiar with R  and basics of statistics


**Description:**

Main focus: "R&Bioconductor and Multiple hypothesis testing"

This tutorial session will provide a brief introduction to Data Science and data manipulation using R and the Jupyter notebooks.

Program:


- R&Bioconductor
- Basic data types, Classes and Methods
- Useful R packages and functions
- Multiple hypothesis Testing (FDR, Bonferroni and Benjamini correction)
- Examples (exercise notebook)

---

## Outline/ content


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


<a href="{{site.url}}{{site.baseurl}}/index.html" class="float" download>
<i class="fa fa-home my-float"></i>
</a>
