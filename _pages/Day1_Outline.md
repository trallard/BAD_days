---
layout: page
permalink: /1-outline/
title: Day 1 outline
category: module
---

# BAD Day 1

**Description**

## Outline
This tutorial session will provide a brief introduction to Data Science and data manipulation using R and the Jupyter notebooks.

Briefly, the topics covered are:

* A brief introduction to Data Science
* R for exploratory data Analysis
* Hands on data: how to produce good graphics

The hands on components of this module are contained in the following tutorials:

1. [Day 1 tutorial]({{site.url}}{{site.baseurl}}/tutorial)
2. [Comments on probability](https://github.com/trallard/BAD_days/blob/master/Day1/Comments.ipynb)


{% assign pages_list = site.pages %}
{% for node in pages_list %}
  {% if node.title != null %}
    {% if node.layout == "default" %}
    {%if node.tags %}


      {% for tag in node.tags %}
      {% if tag == 'Day1' %}
      <!-- Note you need to prepend the site.baseurl always-->
        <a class="sidebar-nav-item{% if page.url == node.url %} active{% endif %}"
        href="{{site.baseurl}}{{ node.url }}">{{ node.title }}</a>
      {% endif %}
      {% endfor %}

    {% endif %}
    {%endif%}
  {% endif %}
{% endfor %}
